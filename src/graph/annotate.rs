use std::collections::{HashMap, HashSet};
use std::convert::TryFrom;
use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};

use anyhow::Result;
use bio::data_structures::interval_tree::ArrayBackedIntervalTree;
use bio::io::gff;
use bio::io::gff::GffType::GFF3;
use clap::Parser;
use flate2::bufread::GzDecoder;
use itertools::Itertools;
use noodles::bcf::header::StringMap;
use noodles::bcf::Reader;
use noodles::vcf::record::info::field::{Key, Value};
use noodles::vcf::{Header, Record};
use noodles::{bcf, vcf};
use serde::Serialize;

use crate::cli::{EventId, GraphStorage};
use crate::graph::EdgeType::Neighbour;
use crate::graph::{CIRCLE_LENGTH_KEY, CIRCLE_SEGMENT_COUNT_KEY, SUPPORT_KEY};

#[derive(Parser)]
pub(crate) struct AnnotateArgs {
    /// Input graph in msgpack format
    #[clap(parse(from_os_str))]
    graph: PathBuf,

    /// VCF/BCF file containing filtered and processed breakend events for the circles described in the graph
    #[clap(parse(from_os_str))]
    breakend_vcf: PathBuf,

    /// Output table with annotated segments in tsv format
    #[clap(parse(from_os_str))]
    output: PathBuf,

    /// GFF3 file with annotations with respect to the reference sequence
    #[clap(parse(from_os_str))]
    annotation: PathBuf,
}

type Annotation = HashMap<String, ArrayBackedIntervalTree<u64, gff::Record>>;
fn read_gff3<P: AsRef<Path> + std::fmt::Debug>(
    path: P,
    filter: fn(&gff::Record) -> bool,
) -> Result<Annotation> {
    let mut annotations = gff::Reader::new(
        File::open(path).map(BufReader::new).map(GzDecoder::new).unwrap(),
        GFF3,
    );
    let mut trees: HashMap<String, ArrayBackedIntervalTree<u64, gff::Record>> = HashMap::new();
    annotations
        .records()
        .flatten()
        .filter(filter)
        .for_each(|r| {
            let tree = trees
                .entry(r.seqname().to_string())
                .or_insert_with(ArrayBackedIntervalTree::new);
            tree.insert(*r.start()..*r.end(), r);
        });
    trees.values_mut().for_each(|tree| tree.index());
    Ok(trees)
}

pub(crate) fn main_annotate(args: AnnotateArgs) -> Result<()> {
    eprintln!("Reading graph");
    let graph = GraphStorage::from_path(&args.graph)?;

    eprintln!("Reading breakend records");
    let (mut reader, header, string_map) = read_bcf(&args)?;

    let event_key: Key = "EVENT".parse()?;
    let event_records = group_event_records(&mut reader, &header, &string_map, &event_key);

    eprintln!("Reading annotation");
    let annotations = read_gff3(&args.annotation, |record| {
        record.feature_type() == "gene" || record.feature_type() == "exon"
    })?;

    //let mut segment_tables: HashMap<EventId, Vec<_>> = HashMap::with_capacity(event_records.len());

    eprintln!("Building circle table");
    let circle_table = graph
        .valid_paths
        .into_iter()
        .flat_map(|(graph_id, circles)| {
            circles
                .into_iter()
                .enumerate()
                .map(|(circle_id, circle)| {
                    for (from, to, edge_info) in circle.edges {
                        // skip non-neighbour edges or those without coverage
                        if !edge_info.edge_type.contains(Neighbour) || edge_info.coverage <= 1e-4 {
                            continue;
                        }
                        // skip edges that cross different contigs
                        if from.0 != to.0 {
                            continue;
                        }
                        let chrom = from.0;
                        if let Some(annot) = annotations.get(&format!("chr{}", chrom + 1)) {
                            let (from, to) = (from.1.min(to.1), from.1.max(to.1));
                            let foo = annot
                                .find(from as u64..to as u64)
                                .iter()
                                .map(|entry| entry.data())
                                .cloned()
                                .collect_vec();
                            dbg!(&foo);
                        }
                    }
                    let event_name = format!("graph_{}_circle_{}", graph_id, circle_id);
                    let records = event_records.get(&event_name).expect("Event not found");
                    get_varlociraptor_info(records, &event_name, circle_id)
                })
                .collect_vec()
        })
        .collect_vec();
    dbg!(&circle_table);
    Ok(())
}

fn get_varlociraptor_info(
    records: &[Record],
    event_name: &str,
    circle_id: usize,
) -> CircleTableInfo {
    let num_exons = 0;
    let genes = vec![];
    let r = &records[0];
    dbg!(&r);
    let info = CircleTableInfo {
        circle_id,
        length: r
            .info()
            .get(&*CIRCLE_LENGTH_KEY)
            .map(|f| match f.value() {
                Value::Integer(i) => usize::try_from(*i).expect("Invalid circle length"),
                _ => panic!("Expected integer value for circle length"),
            })
            .expect("'CircleLength' info field not found"),
        segment_count: r
            .info()
            .get(&*CIRCLE_SEGMENT_COUNT_KEY)
            .map(|f| match f.value() {
                Value::Integer(i) => usize::try_from(*i).expect("Invalid segment count"),
                _ => panic!("Expected integer value for segment count"),
            })
            .expect("'CircleSegmentCount' info field not found"),
        score: r
            .info()
            .get(&*SUPPORT_KEY)
            .map(|f| match f.value() {
                Value::Integer(i) => usize::try_from(*i).expect("Invalid score"),
                _ => panic!("Expected integer value for score"),
            })
            .expect("'Support' info field not found"),
        num_exons,
        genes,
    };
    dbg!(&info);
    info
}

fn read_bcf(args: &AnnotateArgs) -> Result<(Reader<BufReader<File>>, Header, StringMap)> {
    let mut reader = File::open(&args.breakend_vcf)
        .map(BufReader::new)
        .map(bcf::Reader::new)?;
    let _ = reader.read_file_format()?;
    let raw_header = reader.read_header()?;
    let header: vcf::Header = raw_header.parse()?;
    let string_map = raw_header.parse()?;
    Ok((reader, header, string_map))
}

fn group_event_records(
    reader: &mut Reader<BufReader<File>>,
    header: &Header,
    string_map: &StringMap,
    event_key: &Key,
) -> HashMap<String, Vec<Record>> {
    let event_records = reader
        .records()
        .flatten()
        .flat_map(|r| r.try_into_vcf_record(&header, &string_map))
        .map(|r| {
            (
                r.info()
                    .get(event_key)
                    .map(|v| v.value().to_string())
                    .unwrap_or_else(|| "UNKNOWN_EVENT".to_string()),
                r,
            )
        })
        .into_group_map();
    event_records
}

#[derive(Serialize, Debug)]
struct CircleTableInfo {
    circle_id: usize,
    length: usize,
    segment_count: usize,
    score: usize,
    num_exons: usize,
    genes: Vec<String>,
}
