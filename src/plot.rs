use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::path::PathBuf;
use std::str::FromStr;

use anyhow::{Context, Result};
use bam::pileup::AlnType;
use bam::record::{PCR_OR_OPTICAL_DUPLICATE, RECORD_FAILS_QC, RECORD_UNMAPPED, SECONDARY};
use clap::Parser;
use itertools::Itertools;
use plotly::common::{Line, Marker, Mode, Side, Title};
use plotly::layout::{Axis, BarMode, Shape, ShapeLine};
use plotly::{Bar, ImageFormat, Layout, Plot, Rgb, Scatter};
use strum::IntoStaticStr;

use crate::util::split_reads_for_region;

// use crate::time_it;

#[derive(Parser)]
pub(crate) struct PlotArgs {
    /// Input BAM file
    #[clap(parse(from_os_str))]
    input: PathBuf,

    #[clap(long, short)]
    region: Region,

    #[clap(long, short, parse(from_os_str))]
    output: Option<PathBuf>,

    #[clap(short, long, default_value = "0")]
    threads: u16,

    #[clap(short, long)]
    bin_size: Option<usize>,

    #[clap(short, long)]
    flank: Option<u32>,
}

pub(crate) fn main(args: PlotArgs) -> Result<()> {
    let mut reader = bam::IndexedReader::build()
        .additional_threads(args.threads)
        .from_path(args.input)?;
    let header = reader.header();
    let target = header
        .reference_id(&args.region.target)
        .context("No such target name")?;
    let target_end = header.reference_len(target).unwrap_or(u32::MAX);
    let flank = args.flank.unwrap_or(0);
    let original_region = args.region.clone();
    let region = &args.region;
    let plot_region = Region {
        target: region.target.clone(),
        start: region.start.saturating_sub(flank),
        end: target_end.min(region.end.saturating_add(flank)),
    };
    let bam_region = bam::Region::new(target, plot_region.start, plot_region.end);

    let mut bin_size = args
        .bin_size
        .unwrap_or_else(|| (original_region.len() as f64).sqrt().round() as usize);
    let max_bins = 1200;
    if original_region.len() / bin_size > max_bins {
        bin_size = original_region.len() / max_bins;
    }
    let mut detailed_coverage =
        detailed_coverage(&mut reader, &args.region, &bam_region, target_end, bin_size)?;
    let plot = plot(
        &mut detailed_coverage,
        &original_region,
        &plot_region,
        bin_size,
        flank,
    )?;

    if let Some(output) = args.output {
        if let Some("html") = output.extension().and_then(|ext| ext.to_str()) {
            plot.to_html(output);
        } else {
            plot.save(
                output.clone(),
                match output.extension().and_then(|ext| ext.to_str()) {
                    Some("png") => ImageFormat::PNG,
                    Some("svg") => ImageFormat::SVG,
                    Some("pdf") => ImageFormat::PDF,
                    Some(_) => ImageFormat::PNG,
                    None => ImageFormat::PNG,
                },
                1024,
                680,
                1.0,
            );
        }
    } else {
        plot.show()
    }
    Ok(())
}

pub(crate) fn detailed_coverage(
    reader: &mut bam::IndexedReader<File>,
    circle_region: &Region,
    bam_region: &bam::Region,
    ref_len: u32,
    bin_size: usize,
) -> Result<HashMap<Feature, Vec<usize>>> {
    let n_bins = bam_region.len() as usize / bin_size + 1;
    let mut detailed_coverage = HashMap::new();

    let play = 3;
    let split_reads = split_reads_for_region(
        &bam::Region::new(
            bam_region.ref_id(),
            circle_region.start.saturating_sub(play),
            circle_region.end.saturating_add(play).min(ref_len),
        ),
        reader,
    )?;
    let start = circle_region.start as i64;
    let end = circle_region.end as i64;
    enum AtBreakpoint {
        Start,
        End,
        Both,
        False,
    }
    let read_breakpoint = |record: &bam::Record| -> AtBreakpoint {
        let rstart = record.start() as u32;
        let rend = record.calculate_end() as u32;
        let start = (rstart as i64 - start).abs() < play as i64;
        let end = (rend as i64 - end).abs() < play as i64;
        match (start, end) {
            (true, true) => AtBreakpoint::Both,
            (true, false) => AtBreakpoint::Start,
            (false, true) => AtBreakpoint::End,
            _ => AtBreakpoint::False,
        }
    };
    let has_breakpoint_partner = |record: &bam::Record| -> bool {
        let name = std::str::from_utf8(record.name()).unwrap().to_owned();
        if let Some(partners) = split_reads.get(&name) {
            let breakpoint = read_breakpoint(record);
            match breakpoint {
                AtBreakpoint::Start => partners
                    .iter()
                    .map(read_breakpoint)
                    .any(|x| matches!(x, AtBreakpoint::End)),
                AtBreakpoint::End => partners
                    .iter()
                    .map(read_breakpoint)
                    .any(|x| matches!(x, AtBreakpoint::Start)),
                AtBreakpoint::Both => true,
                AtBreakpoint::False => false,
            }
        } else {
            false
        }
    };
    let inside: HashSet<_> = split_reads
        .iter()
        .filter(|(_, records)| records.len() > 1)
        .map(|(name, _)| name)
        .collect();

    let filter = |record: &bam::Record| -> bool {
        record
            .flag()
            .no_bits(RECORD_FAILS_QC | PCR_OR_OPTICAL_DUPLICATE | RECORD_UNMAPPED | SECONDARY)
    };

    let mut reads = reader.fetch_by(bam_region, filter)?;
    let pileups = bam::Pileup::new(&mut reads);

    let mut read_depths = vec![0usize; n_bins];
    for pileup in pileups {
        let pileup = pileup?;
        let ref_pos = pileup.ref_pos();
        if ref_pos >= bam_region.end() {
            break;
        }
        let idx = (ref_pos.saturating_sub(bam_region.start()) as usize) / bin_size;
        read_depths[idx] = read_depths[idx].saturating_add(pileup.entries().len());
        for entry in pileup.entries() {
            // AlnType doesn't implement Hash, so we have our own Feature enum
            let aln = match entry.aln_type() {
                AlnType::Deletion => Feature::Deletion,
                AlnType::Match => Feature::Match,
                AlnType::Insertion(_) => Feature::Insertion,
            };

            let record = entry.record();
            if record.flag().is_supplementary() || record.tags().get(b"SA").is_some() {
                let name = std::str::from_utf8(record.name())?.to_owned();
                detailed_coverage
                    .entry(match has_breakpoint_partner(record) {
                        true => Feature::SplitReadAtBreakpoints,
                        false => {
                            if inside.contains(&name) {
                                Feature::SplitReadIn
                            } else {
                                Feature::SplitReadOut
                            }
                        }
                    })
                    .or_insert_with(|| vec![0usize; n_bins])[idx] += 1;
            } else {
                detailed_coverage
                    .entry(aln)
                    .or_insert_with(|| vec![0usize; n_bins])[idx] += 1;
            }
        }
    }
    detailed_coverage.insert(Feature::Any, read_depths);
    Ok(detailed_coverage)
}

pub(crate) fn plot(
    coverage: &mut HashMap<Feature, Vec<usize>>,
    target_region: &Region,
    plot_region: &Region,
    bin_size: usize,
    _flank: u32,
) -> Result<Plot> {
    let mut plot = Plot::new();
    let read_depths = coverage.remove(&Feature::Any).unwrap();
    let max_read_depth = *read_depths.iter().max().unwrap_or(&0);

    let xx: Vec<u32> = (plot_region.start..plot_region.end)
        .step_by(bin_size)
        .collect();

    for (feature, counts) in coverage
        .iter()
        .sorted_by_key(|&(&feature, _counts)| feature)
    {
        let bar = Bar::new(
            xx.clone(),
            counts
                .iter()
                .zip(read_depths.iter())
                .map(|(&count, &sum)| count as f64 / sum as f64)
                .collect_vec(),
        )
        .name(feature.into())
        .width((bin_size as f64 * 0.95).round() as usize)
        // from https://personal.sron.nl/~pault/
        .marker(Marker::new().color(match feature {
            Feature::Match => Rgb::new(0, 119, 187),         // "blue"
            Feature::Deletion => Rgb::new(0, 153, 136),      // "teal"
            Feature::Insertion => Rgb::new(51, 187, 238),    // "cyan"
            Feature::SplitReadIn => Rgb::new(238, 119, 51),  // "orange"
            Feature::SplitReadOut => Rgb::new(238, 51, 119), // "magenta"
            Feature::SplitReadAtBreakpoints => Rgb::new(204, 51, 17), // "red"
            Feature::Any => Rgb::new(0, 0, 0),
        }))
        //.opacity(0.9)
        .x_axis("x")
        .y_axis("y");
        plot.add_trace(bar);
    }

    let y = read_depths;
    let read_depth_trace = Scatter::new(xx, y)
        .x_axis("x")
        .y_axis("y2")
        .marker(Marker::new().color(Rgb::new(187, 187, 187)).size(6)) // grey
        .line(Line::new().width(4.))
        .name("read depth")
        .mode(Mode::LinesMarkers);
    plot.add_trace(read_depth_trace);
    let layout = Layout::new()
        .template("plotly_white")
        .bar_mode(BarMode::Stack)
        .x_axis(
            Axis::new()
                .title(Title::new("target position"))
                .range(vec![plot_region.start as f64, plot_region.end as f64]),
        )
        .y_axis(
            Axis::new()
                .title(Title::new("ratio of bases"))
                .anchor("x")
                .range(vec![0., 1.])
                .side(Side::Left),
        )
        .y_axis2(
            Axis::new()
                .title(Title::new("read depth"))
                .overlaying("y")
                .anchor("x")
                .range(vec![0., max_read_depth as f64])
                .side(Side::Right),
        )
        .title(Title::new(&format!(
            "Coverage for region {}:{}-{}, length: {} (bin_size: {})",
            &target_region.target,
            target_region.start,
            target_region.end,
            target_region.end - target_region.start,
            bin_size
        )))
        .shapes(vec![
            Shape::new()
                .line(ShapeLine::new().dash("dot"))
                .x0(target_region.start as f64 - bin_size as f64 / 2.)
                .x1(target_region.start as f64 - bin_size as f64 / 2.)
                .y0(0.)
                .y1(1.),
            Shape::new()
                .line(ShapeLine::new().dash("dot"))
                .x0(target_region.end as f64 - bin_size as f64 / 2.)
                .x1(target_region.end as f64 - bin_size as f64 / 2.)
                .y0(0.)
                .y1(1.),
        ]);
    plot.set_layout(layout);
    Ok(plot)
}
//
// /// Extract split reads from BAM.
// /// Here, split read is defined as: record that is either supplementary or has the `SA` tag set,
// /// but is not secondary.
// fn split_reads(reads: &mut RegionViewer<File>) -> Result<HashMap<String, Vec<bam::Record>>> {
//     let mut reads = reads
//         .into_iter()
//         .flatten()
//         .filter(|read| {
//             let flag = read.flag();
//             let tags = read.tags();
//             flag.no_bits(SECONDARY | PCR_OR_OPTICAL_DUPLICATE | RECORD_FAILS_QC | RECORD_UNMAPPED)
//                 && (flag.is_supplementary() || tags.get(b"SA").is_some())
//         })
//         .filter_map(|read| {
//             let name = std::str::from_utf8(read.name()).map(|name| name.to_owned());
//             name.map(|name| (name, read)).ok()
//         })
//         .into_group_map();
//
//     reads = reads
//         .into_iter()
//         .filter(|(_, values)| !values.is_empty())
//         .collect();
//     Ok(reads)
// }

#[derive(IntoStaticStr, Debug, Clone, Hash, Copy, Eq, PartialEq, Ord, PartialOrd)]
pub(crate) enum Feature {
    Any,
    Match,
    Deletion,
    Insertion,
    SplitReadIn,
    SplitReadOut,
    SplitReadAtBreakpoints,
}

#[derive(Debug, Clone)]
pub struct Region {
    pub(crate) target: String,
    pub(crate) start: u32,
    pub(crate) end: u32,
}

impl Region {
    pub(crate) fn new(target: String, start: u32, end: u32) -> Self {
        Self { target, start, end }
    }

    pub(crate) fn len(&self) -> usize {
        (self.end - self.start) as usize
    }
}

impl FromStr for Region {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let (target, range) = s.split_once(':').context("No ':' in region string")?;
        let (start, end) = range.split_once('-').context("No '-' in region string")?;
        let start = start.parse::<u32>()?;
        let end = end.parse::<u32>()?;
        Ok(Region {
            target: target.into(),
            start,
            end,
        })
    }
}
