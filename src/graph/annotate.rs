use std::collections::HashMap;
use std::convert::TryFrom;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::Result;
use bio::data_structures::interval_tree::ArrayBackedIntervalTree;
use bio::io::gff;
use bio::io::gff::GffType::GFF3;
use clap::Parser;
use flate2::bufread::MultiGzDecoder;
use itertools::Itertools;
use noodles::bcf::header::StringMap;
use noodles::bcf::Reader;
use noodles::vcf::record::info::field::{Key, Value};
use noodles::vcf::{Header, Record};
use noodles::{bcf, vcf};
use serde::Serialize;

use crate::cli::{CircleId, EventId, GraphStorage};
use crate::common::ReferenceId;
use crate::graph::EdgeType::{Neighbour, Split};
use crate::graph::{Cycle, Position};

#[derive(Parser)]
pub(crate) struct AnnotateArgs {
    /// Input graph in msgpack format
    #[clap(parse(from_os_str))]
    graph: PathBuf,

    /// VCF/BCF file containing filtered and processed breakend events for the circles described in the graph
    #[clap(parse(from_os_str))]
    breakend_vcf: PathBuf,

    /// (b)gzipped GFF3 file containing annotations with respect to the reference sequence
    #[clap(long, parse(from_os_str))]
    annotation: PathBuf,

    /// Path for the overview circle table
    #[clap(long, parse(from_os_str))]
    circle_table: PathBuf,

    /// Output directory for detailed per-segment information for each circle
    #[clap(long, parse(from_os_str))]
    segment_tables: PathBuf,
}

type Annotation = HashMap<String, ArrayBackedIntervalTree<u64, gff::Record>>;
fn read_gff3<P: AsRef<Path> + std::fmt::Debug>(
    path: P,
    filter: fn(&gff::Record) -> bool,
) -> Result<Annotation> {
    let mut annotations = gff::Reader::new(
        File::open(path)
            .map(BufReader::new)
            .map(MultiGzDecoder::new)
            .unwrap(),
        GFF3,
    );
    let mut trees: HashMap<String, ArrayBackedIntervalTree<u64, gff::Record>> = HashMap::new();
    annotations
        .records()
        .map(|r| r.unwrap())
        .filter(filter)
        .for_each(|r| {
            let tree = trees
                .entry(r.seqname().to_string())
                .or_insert_with(ArrayBackedIntervalTree::new);
            let (start, end) = (r.start().min(r.end()), r.start().max(r.end()));
            tree.insert(*start..*end, r);
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

    std::fs::create_dir_all(&args.segment_tables)?;

    eprintln!("Building circle table");
    let circle_table = graph
        .valid_paths
        .into_iter()
        .flat_map(|(graph_id, circles)| {
            circles
                .into_iter()
                .enumerate()
                .filter_map(|(circle_id, circle)| {
                    let event_name = format!("graph_{}_circle_{}", graph_id, circle_id);
                    if let Some(records) = event_records.get(&event_name) {
                        let chrom = format!("{}", records[0].chromosome());
                        let chrom = chrom.strip_prefix("chr").unwrap_or(&chrom);
                        let segment_annotations = segment_annotation(&annotations, circle, chrom);
                        let varlociraptor_annotations = varlociraptor_info(records);

                        let mut segment_writer =
                            csv::WriterBuilder::new().delimiter(b'\t').from_writer(
                                File::create(
                                    &args
                                        .segment_tables
                                        .join(format!("{}_segments.tsv", event_name)),
                                )
                                .map(BufWriter::new)
                                .unwrap(),
                            );
                        let entry = (
                            graph_id,
                            circle_id,
                            varlociraptor_annotations,
                            segment_annotations,
                        );
                        write_segment_table_into(&mut segment_writer, &entry).unwrap();
                        Some(entry)
                    } else {
                        eprintln!("WARNING: Event {} not found in breakend VCF", event_name);
                        None
                    }
                })
                .collect_vec()
        })
        .collect_vec();

    write_circle_table(
        File::create(&args.circle_table).map(BufWriter::new)?,
        &circle_table,
    )?;
    Ok(())
}

type Segment = (ReferenceId, Position, ReferenceId, Position);
type CircleTableEntry = (
    EventId,
    CircleId,
    CircleInfo,
    Vec<(Segment, Option<SegmentAnnotation>)>,
);

fn write_circle_table<W: Write>(writer: W, circle_table: &[CircleTableEntry]) -> Result<()> {
    let mut writer: csv::Writer<_> = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(writer);
    for (graph_id, circle_id, circle_info, segment_annotations) in circle_table {
        let (num_exons, gene_ids, gene_names, num_split_reads) = segment_annotations.iter().fold(
            (0, vec![], vec![], 0),
            |(mut n_ex, mut gids, mut gnames, mut n_splits), (_, sa)| {
                if let Some(sa) = sa {
                    match sa {
                        SegmentAnnotation::Segment {
                            num_exons,
                            gene_ids,
                            gene_names,
                            num_split_reads,
                            coverage: _,
                        } => {
                            n_ex += num_exons;
                            gids.extend(gene_ids);
                            gnames.extend(gene_names);
                            n_splits += if circle_info.segment_count == 1 {
                                // if there's only one segment, the split reads loop from one end of the segment to its other end
                                *num_split_reads
                            } else {
                                // otherwise, split reads are described by junction edges and counted below
                                0
                            };
                        }
                        SegmentAnnotation::Junction { num_split_reads } => {
                            n_splits += num_split_reads;
                        }
                    }
                }
                (n_ex, gids, gnames, n_splits)
            },
        );
        let gene_ids = gene_ids.into_iter().unique().join(",");
        let gene_names = gene_names.into_iter().unique().join(",");
        let cti = FlatCircleTableInfo {
            graph_id: *graph_id,
            circle_id: *circle_id,
            circle_length: circle_info.length,
            segment_count: circle_info.segment_count,
            score: circle_info.score,
            num_split_reads,
            num_exons,
            gene_ids,
            gene_names,
        };
        writer.serialize(cti)?;
    }
    Ok(())
}

fn write_segment_table_into<W: Write>(
    writer: &mut csv::Writer<W>,
    circle_entry: &CircleTableEntry,
) -> Result<()> {
    let (graph_id, circle_id, circle_info, segment_annotations) = circle_entry;
    #[derive(Serialize, Debug)]
    struct SegmentEntry {
        graph_id: usize,
        circle_id: usize,
        circle_length: usize,
        segment_count: usize,
        score: usize,
        kind: Option<String>,
        target_from: String,
        from: Position,
        target_to: String,
        to: Position,
        length: Option<Position>,
        num_exons: Option<usize>,
        gene_ids: Option<String>,
        gene_names: Option<String>,
        coverage: Option<f64>,
        num_split_reads: Option<usize>,
    }
    for ((from_ref_id, from, to_ref_id, to), annotation) in segment_annotations {
        let (segment_type, num_exons, gene_names, gene_ids, num_split_reads, coverage) =
            if let Some(annotation) = annotation.as_ref() {
                match annotation {
                    SegmentAnnotation::Segment {
                        num_exons,
                        gene_ids,
                        gene_names,
                        num_split_reads,
                        coverage,
                    } => (
                        Some("coverage".to_string()),
                        Some(*num_exons),
                        Some(gene_names.join(",")),
                        Some(gene_ids.join(",")),
                        Some(*num_split_reads),
                        Some(*coverage),
                    ),
                    SegmentAnnotation::Junction { num_split_reads } => (
                        Some("split".to_string()),
                        None,
                        None,
                        None,
                        Some(*num_split_reads),
                        None,
                    ),
                }
            } else {
                (None, None, None, None, None, None)
            };
        let is_coverage_segment = if let Some(k) = segment_type.as_ref() {
            k == "coverage"
        } else {
            false
        };
        let entry = SegmentEntry {
            graph_id: *graph_id,
            circle_id: *circle_id,
            circle_length: circle_info.length,
            segment_count: circle_info.segment_count,
            score: circle_info.score,
            kind: segment_type,
            target_from: format!("todo: {}", from_ref_id),
            from: *from,
            target_to: format!("todo: {}", to_ref_id),
            to: *to,
            length: if is_coverage_segment {
                Some(from.abs_diff(*to))
            } else {
                None
            },
            num_exons,
            gene_ids,
            gene_names,
            coverage,
            num_split_reads,
        };
        writer.serialize(entry)?;
    }
    Ok(())
}

fn segment_annotation(
    annotations: &Annotation,
    circle: Cycle,
    chrom: &str,
) -> Vec<(Segment, Option<SegmentAnnotation>)> {
    circle
        .edges
        .iter()
        .map(|&((from_ref_id, from), (to_ref_id, to), edge_info)| {
            let is_coverage_segment = edge_info.edge_type.contains(Neighbour)
                && edge_info.coverage > 1e-4
                && from_ref_id == to_ref_id;
            if is_coverage_segment {
                if let Some(annot) = annotations
                    .get(chrom)
                    .or_else(|| annotations.get(&format!("chr{}", chrom)))
                {
                    let (from, to) = (from.min(to), from.max(to));
                    let (num_exons, gene_ids, gene_names) = annot
                        .find(from as u64..to as u64)
                        .iter()
                        .map(|entry| entry.data())
                        .filter_map(|data| match data.feature_type() {
                            "gene" => Some((
                                0,
                                data.attributes().get("gene_id").cloned(),
                                data.attributes().get("gene_name").cloned(),
                            )),
                            "exon" => Some((
                                1,
                                data.attributes().get("gene_id").cloned(),
                                data.attributes().get("gene_name").cloned(),
                            )),
                            _ => None,
                        })
                        .fold(
                            (0, vec![], vec![]),
                            |(mut num_exons, mut gene_ids, mut gene_names),
                             (n_exons, gene_id, gene_name)| {
                                num_exons += n_exons;
                                gene_ids.extend(gene_id);
                                gene_names.extend(gene_name);
                                (num_exons, gene_ids, gene_names)
                            },
                        );
                    (
                        (from_ref_id, from, to_ref_id, to),
                        Some(SegmentAnnotation::Segment {
                            num_exons,
                            gene_ids,
                            gene_names,
                            coverage: edge_info.coverage,
                            num_split_reads: edge_info.num_split_reads,
                        }),
                    )
                } else {
                    ((from_ref_id, from, to_ref_id, to), None)
                }
            } else {
                // not a coverage segment but a split junction
                if edge_info.edge_type.contains(Split) {
                    (
                        (from_ref_id, from, to_ref_id, to),
                        Some(SegmentAnnotation::Junction {
                            num_split_reads: edge_info.num_split_reads,
                        }),
                    )
                } else {
                    ((from_ref_id, from, to_ref_id, to), None)
                }
            }
        })
        .collect_vec()
}

fn varlociraptor_info(records: &[Record]) -> CircleInfo {
    let r = &records[0];
    let circle_length_key = "CircleLength".parse::<Key>().unwrap();
    let circle_segment_count_key = "CircleSegmentCount".parse::<Key>().unwrap();
    let support_key = "Support".parse::<Key>().unwrap();

    let length = r
        .info()
        .get(&circle_length_key)
        .map(|f| match f.value() {
            Value::Integer(i) => usize::try_from(*i).expect("Invalid circle length"),
            _ => panic!("Expected integer value for circle length"),
        })
        .expect("'CircleLength' info field not found");
    let segment_count = r
        .info()
        .get(&circle_segment_count_key)
        .map(|f| match f.value() {
            Value::Integer(i) => usize::try_from(*i).expect("Invalid segment count"),
            _ => panic!("Expected integer value for segment count"),
        })
        .expect("'CircleSegmentCount' info field not found");
    let score = r
        .info()
        .get(&support_key)
        .map(|f| match f.value() {
            Value::Integer(i) => usize::try_from(*i).expect("Invalid score"),
            _ => panic!("Expected integer value for score"),
        })
        .expect("'Support' info field not found");

    CircleInfo {
        length,
        segment_count,
        score,
    }
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
        // .map(|r| {
        //     r.try_into_vcf_record(header, string_map)
        //         .unwrap_or_else(|_| panic!("{:?}:{:?}", &r.chromosome_id(), &r.position()))
        // })
        .flat_map(|r| r.try_into_vcf_record(header, string_map))
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

#[derive(Serialize, Debug, Clone, Copy)]
struct CircleInfo {
    length: usize,
    segment_count: usize,
    score: usize,
}

#[derive(Debug)]
enum SegmentAnnotation {
    Segment {
        num_exons: usize,
        gene_ids: Vec<String>,
        gene_names: Vec<String>,
        num_split_reads: usize,
        coverage: f64,
    },
    Junction {
        num_split_reads: usize,
    },
}

#[derive(Serialize, Debug, Default)]
struct CircleAnnotation {
    num_exons: usize,
    gene_ids: String,
    gene_names: String,
}

// serialization of this with csv doesn't work when the header is derived
// see https://github.com/BurntSushi/rust-csv/issues/188 and https://github.com/BurntSushi/rust-csv/pull/197
// #[derive(Serialize, Debug)]
// struct CircleTableInfo {
//     graph_id: usize,
//     circle_id: usize,
//     circle_info: CircleInfo,
//     annotation_info: CircleAnnotation,
// }
#[derive(Serialize, Debug)]
struct FlatCircleTableInfo {
    graph_id: usize,
    circle_id: usize,
    circle_length: usize,
    segment_count: usize,
    score: usize,
    num_exons: usize,
    gene_ids: String,
    gene_names: String,
    num_split_reads: usize,
}
