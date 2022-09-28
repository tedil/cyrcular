use std::cmp::Ordering;
use std::collections::HashMap;
use std::hash::{Hash, Hasher};
use std::io::{Read, Seek};

use anyhow::Result;
use bam::record::{PCR_OR_OPTICAL_DUPLICATE, RECORD_FAILS_QC, RECORD_UNMAPPED, SECONDARY};
use bam::{Header, Region};
use itertools::Itertools;
use noodles::vcf::header::record::value::{map::Contig, Map};

use crate::common::SplitReadStorage;
use crate::graph::{CIRCLE_LENGTH_KEY, CIRCLE_SEGMENT_COUNT_KEY, NUM_SPLIT_READS_KEY, SUPPORT_KEY};

// bam::Region doesn't implement Eq or Hash, so we'll have to wrap that for now
#[derive(Clone, Debug)]
pub struct HRegion(pub Region);

impl Hash for HRegion {
    fn hash<H: Hasher>(&self, state: &mut H) {
        state.write_u32(self.0.ref_id());
        state.write_u32(self.0.start());
        state.write_u32(self.0.end());
        state.finish();
    }
}

impl PartialEq for HRegion {
    fn eq(&self, other: &Self) -> bool {
        self.0.ref_id() == other.0.ref_id()
            && self.0.start() == other.0.start()
            && self.0.end() == other.0.end()
    }
}

impl Eq for HRegion {}

impl PartialOrd for HRegion {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(
            self.0
                .ref_id()
                .cmp(&other.0.ref_id())
                .then(self.0.start().cmp(&other.0.start()))
                .then(self.0.end().cmp(&other.0.end())),
        )
    }
}

impl Ord for HRegion {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other).unwrap()
    }
}

/// Measure time to execute some code fragment.
///
/// Inspired by
///   <https://notes.iveselov.info/programming/time_it-a-case-study-in-rust-macros>
#[macro_export]
macro_rules! time_it {
    ($context:literal, $s:expr) => {{
        use std::io::Write;
        eprint!("{:25}: ", $context);
        std::io::stderr().flush().unwrap();
        let timer = std::time::Instant::now();
        let r = $s;
        eprintln!("{:08.3?}", timer.elapsed());
        r
    }};
    ($context:expr, $s:expr) => {{
        use std::io::Write;
        eprint!("{:25}: ", $context);
        std::io::stderr().flush().unwrap();
        let timer = std::time::Instant::now();
        let r = $s;
        eprintln!("{:08.3?}", timer.elapsed());
        r
    }};
}

#[inline]
pub(crate) fn default_filter(record: &bam::Record) -> bool {
    record
        .flag()
        .no_bits(RECORD_FAILS_QC | PCR_OR_OPTICAL_DUPLICATE | RECORD_UNMAPPED | SECONDARY)
}

#[inline]
pub(crate) fn is_split_read(record: &bam::Record) -> bool {
    record.flag().is_supplementary() || record.tags().get(b"SA").is_some()
}

#[inline]
pub(crate) fn split_read_filter(record: &bam::Record) -> bool {
    is_split_read(record) && default_filter(record)
}

#[allow(dead_code)]
pub(crate) fn split_reads_for_region(
    region: &bam::Region,
    reader: &mut bam::IndexedReader<std::fs::File>,
) -> Result<HashMap<String, Vec<bam::Record>>> {
    let iter = reader.fetch_by(region, split_read_filter)?;
    Ok(iter
        .flatten()
        .map(|record| {
            (
                std::str::from_utf8(record.name()).unwrap().to_owned(),
                record.clone(),
            )
        })
        .into_group_map())
}

pub(crate) fn split_reads<R: Seek + Read>(
    reader: &mut bam::IndexedReader<R>,
) -> Result<SplitReadStorage> {
    let iter = reader.full_by(split_read_filter);
    Ok(iter
        .flatten()
        .map(|record| {
            (
                std::str::from_utf8(record.name()).unwrap().to_owned(),
                record.clone(),
            )
        })
        .into_group_map())
}

pub fn vcf_header_from_bam(bam_header: &Header) -> noodles::vcf::Header {
    let builder = bam_header
        .reference_names()
        .iter()
        .zip(bam_header.reference_lengths())
        .fold(
            noodles::vcf::Header::builder(),
            |builder, (name, length)| {
                let mut contig = Map::<Contig>::new(name.parse().unwrap());
                *contig.length_mut() = Some(*length as usize);
                builder.add_contig(contig)
            },
        );
    use noodles::vcf::header::{
        info::Key,
        record::value::{map::Filter, map::Info},
    };
    let builder = builder
        .add_info(Map::<Info>::from(Key::SvType))
        .add_info(Map::<Info>::from(Key::MateBreakendIds))
        .add_info(Map::<Info>::from(Key::BreakendEventId))
        .add_info(Map::<Info>::from((*SUPPORT_KEY).clone()))
        .add_info(Map::<Info>::from((*CIRCLE_LENGTH_KEY).clone()))
        .add_info(Map::<Info>::from((*CIRCLE_SEGMENT_COUNT_KEY).clone()))
        .add_info(Map::<Info>::from((*NUM_SPLIT_READS_KEY).clone()))
        .add_filter(Map::<Filter>::pass());
    builder.build()
}
