use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::Result;
use clap::Parser;
use itertools::Itertools;
use noodles::bcf::header::StringMaps;
use noodles::bcf::Reader;
use noodles::bgzf::reader::Reader as BgzfReader;
use noodles::vcf::header::info::Key;
use noodles::vcf::record::info::field::Value;
use noodles::vcf::{Header, Record};
use noodles::{bcf, vcf};
use ordered_float::OrderedFloat;
use serde::{Deserialize, Serialize};
use struct_field_names_as_array::FieldNamesAsArray;

use crate::cli::CircleId;
use crate::graph::annotate::{
    event_id, AnnotatedCircle, AnnotatedGraph, GraphId, Reference, SegmentAnnotation,
};
use crate::graph::Position;

#[derive(Parser)]
pub(crate) struct TableArgs {
    /// Path to annotated graph
    #[arg()]
    annotated_graph: PathBuf,

    /// Path to breakend records in BCF format
    #[arg()]
    records: PathBuf,

    /// Path to reference
    #[arg(long)]
    reference: PathBuf,
    
    /// Comma-separated list of (lowercase) event names that have a PROB_<event_name> INFO field in the BCF.
    /// Usually, these will be events that have been specified in the varlociraptor calling scenario YAML file.
    #[arg(long, value_delimiter = ',', num_args = 1..)]
    event_names: Vec<String>,

    /// Path for the overview circle table output
    #[arg(long)]
    circle_table: PathBuf,

    /// Path for the output directory for the detailed per-segment information for each circle
    #[arg(long)]
    segment_tables: PathBuf,
}

pub(crate) fn main_table(args: TableArgs) -> Result<()> {
    let annotated_graph = AnnotatedGraph::from_path(&args.annotated_graph)?;
    eprintln!("Reading breakend records");
    let (mut reader, header, string_maps) = read_bcf(&args.records)?;
    let event_records = group_event_records(&mut reader, &header, &string_maps);
    let reference = Reference::from_path(&args.reference)?;
    let prob_info_fields: Vec<String> = event_names_to_prob_names(&args.event_names);
    let (circle_summaries, details) =
        tables_from_annotated_graph(&annotated_graph, &event_records, &reference, &prob_info_fields)?;

    std::fs::create_dir_all(&args.segment_tables)?;
    write_circle_table(
        File::create(&args.circle_table).map(BufWriter::new)?,
        &circle_summaries,
    )?;

    for (event_name, (_, annotated_circle)) in details {
        let mut segment_writer = segment_table_writer(&args, &event_name)
            .unwrap_or_else(|_| panic!("Failed creating file for {}", &event_name));
        write_segment_table_into(&mut segment_writer, &annotated_circle, &reference)?;
    }

    Ok(())
}

pub(crate) type BcfReader = Reader<BgzfReader<BufReader<File>>>;
pub(crate) fn read_bcf<P: AsRef<Path>>(path: P) -> Result<(BcfReader, Header, StringMaps)> {
    let mut reader = File::open(path.as_ref())
        .map(BufReader::new)
        .map(bcf::Reader::new)?;
    let _ = reader.read_file_format()?;
    let raw_header = reader.read_header()?;
    let header: Header = raw_header.parse()?;
    let string_maps: StringMaps = raw_header.parse()?;
    Ok((reader, header, string_maps))
}
pub(crate) fn group_event_records(
    reader: &mut BcfReader,
    header: &Header,
    string_maps: &StringMaps,
) -> HashMap<String, Vec<Record>> {
    let event_key: Key = "EVENT".parse().unwrap(); // guaranteed to exist
    let event_records = reader
        .records()
        .map(|r| r.unwrap())
        .map(|r| r.try_into_vcf_record(header, string_maps).unwrap())
        .map(|r| {
            (
                r.info()
                    .get(&event_key)
                    .map(|v| v.value().expect("empty EVENT?").to_string())
                    .unwrap_or_else(|| "UNKNOWN_EVENT".to_string()),
                r,
            )
        })
        .into_group_map();
    event_records
}

fn event_names_to_prob_names(event_names: &Vec<String>) -> Vec<String> {
    event_names.into_iter()
        .map(|event_name| event_name.to_uppercase())
        .map(|mut e_n| {
            e_n.insert_str(0, "PROB_");
            e_n
        })
        .collect()
}

fn get_value_from_prob_name(prob_name: &String, record: &Record) -> f32 {
    let key = prob_name.parse::<Key>().unwrap();
    record.info()
        .get(&key)
        .map(|f| match f.value() {
            Some(Value::Float(f)) => *f,
            _ => panic!("Expected float value for {}", prob_name),
        })
        .unwrap_or_else(|| panic!("'{}' info field not found", prob_name))
}

fn varlociraptor_info(records: &[Record], prob_info_fields: &Vec<String>) -> VarlociraptorInfo {
    let r = &records[0];

    let prob_absent = get_value_from_prob_name(&String::from("PROB_ABSENT"), r);
    let prob_artifact = get_value_from_prob_name(&String::from("PROB_ARTIFACT"), r);
    let event_probs: HashMap<String, f32> = prob_info_fields
    .into_iter()
    .map(|field_name| {
        (field_name.clone(), get_value_from_prob_name(field_name, r))
    }).collect();
    let af_key = "AF".parse::<vcf::header::format::Key>().unwrap();
    let allele_frequencies = r
        .genotypes()
        .iter()
        .map(|f| f[&af_key].value())
        .collect_vec();
    VarlociraptorInfo {
        prob_absent,
        prob_artifact,
        event_probs,
        af_nanopore: allele_frequencies[0].map(|v| match v {
            vcf::record::genotype::field::Value::Float(f) => *f,
            _ => panic!("Expected float value for AF"),
        }),
    }
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub(crate) struct VarlociraptorInfo {
    prob_absent: f32,
    prob_artifact: f32,
    event_probs: HashMap<String, f32>,
    af_nanopore: Option<f32>,
}

fn tables_from_annotated_graph(
    graph: &AnnotatedGraph,
    event_records: &HashMap<String, Vec<Record>>,
    reference: &Reference,
    prob_info_fields: &Vec<String>,
) -> Result<(
    Vec<FlatCircleTableInfo>,
    HashMap<String, (VarlociraptorInfo, AnnotatedCircle)>,
)> {
    eprintln!("Building circle table");
    let event_details: HashMap<String, _> = graph
        .annotated_circles
        .iter()
        .cloned()
        .flat_map(|(graph_id, circles)| {
            circles
                .into_iter()
                .enumerate()
                .filter_map(|(circle_id, annotated_circle)| {
                    let event_name = format!("graph_{}_circle_{}", graph_id, circle_id);
                    if let Some(records) = event_records.get(&event_name) {
                        Some((event_name, (varlociraptor_info(records, prob_info_fields), annotated_circle)))
                    } else {
                        eprintln!("WARNING: Event {} not found in breakend VCF", &event_name);
                        None
                    }
                })
                .collect_vec()
        })
        .collect();

    let summaries = graph.summaries(reference);
    let circle_table = summaries
        .into_iter()
        .filter_map(|summary| {
            let gene_ids = summary.gene_ids.into_iter().unique().join(",");
            let gene_names = summary.gene_names.into_iter().unique().join(",");
            let regulatory_features = summary.regulatory_features.into_iter().join(",");
            let repeats = summary.repeats.into_iter().join(",");
            let regions = if summary.segment_count == 1 {
                summary.regions.first().cloned().unwrap_or_default()
            } else {
                summary.regions.join(",")
            };
            let event_name = format!("graph_{}_circle_{}", summary.graph_id, summary.circle_id);
            event_details
                .get(&event_name)
                .map(|(varlociraptor_info, _)| FlatCircleTableInfo {
                    event_id: event_id(summary.graph_id, summary.circle_id),
                    graph_id: summary.graph_id,
                    circle_id: summary.circle_id,
                    circle_length: summary.circle_length,
                    segment_count: summary.segment_count,
                    regions,
                    num_split_reads: summary.num_split_reads,
                    num_exons: summary.num_exons,
                    gene_ids,
                    gene_names,
                    regulatory_features,
                    repeats,
                    varlociraptor_event_prob: varlociraptor_info.event_probs.values().sum(),
                    prob_absent: varlociraptor_info.prob_absent,
                    prob_artifact: varlociraptor_info.prob_absent,
                    af_nanopore: varlociraptor_info.af_nanopore,
                })
        })
        .sorted_unstable_by_key(|table_entry| {
            (
                OrderedFloat(1. - table_entry.varlociraptor_event_prob),
                -(table_entry.num_exons as i64),
                OrderedFloat(table_entry.prob_absent),
                table_entry.graph_id,
                table_entry.circle_id,
            )
        })
        .collect_vec();
    Ok((circle_table, event_details))
}

fn segment_table_writer(
    args: &TableArgs,
    event_name: &String,
) -> Result<csv::Writer<BufWriter<File>>> {
    Ok(csv::WriterBuilder::new().delimiter(b'\t').from_writer(
        File::create(
            args.segment_tables
                .join(format!("{}_segments.tsv", event_name)),
        )
        .map(BufWriter::new)?,
    ))
}

fn write_circle_table<W: Write>(writer: W, circle_table: &[FlatCircleTableInfo]) -> Result<()> {
    let mut writer: csv::Writer<_> = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_writer(writer);
    writer.write_record(FlatCircleTableInfo::FIELD_NAMES_AS_ARRAY)?;
    for entry in circle_table {
        writer.serialize(entry)?;
    }
    Ok(())
}

fn write_segment_table_into<W: Write>(
    writer: &mut csv::Writer<W>,
    annotated_circle: &AnnotatedCircle,
    reference: &Reference,
) -> Result<()> {
    #[derive(Serialize, Debug)]
    struct SegmentEntry {
        graph_id: GraphId,
        circle_id: CircleId,
        kind: Option<String>,
        target_from: String,
        from: Position,
        target_to: String,
        to: Position,
        length: Option<Position>,
        num_exons: Option<usize>,
        gene_ids: Option<String>,
        gene_names: Option<String>,
        regulatory_features: Option<String>,
        repeats: Option<String>,
        coverage: Option<f64>,
        num_split_reads: Option<usize>,
        breakpoint_sequence: Option<String>,
    }
    for ((from_ref_id, from, to_ref_id, to), annotation) in annotated_circle.segments() {
        let (
            segment_type,
            num_exons,
            gene_names,
            gene_ids,
            regulatory_features,
            repeats,
            num_split_reads,
            coverage,
            breakpoint_sequence,
        ) = if let Some(annotation) = annotation.as_ref() {
            match annotation {
                SegmentAnnotation::Segment {
                    num_exons,
                    gene_ids,
                    gene_names,
                    regulatory_features,
                    repeats,
                    num_split_reads,
                    coverage,
                    breakpoint_sequence,
                } => (
                    Some("coverage".to_string()),
                    Some(*num_exons),
                    Some(gene_names.join(",")),
                    Some(gene_ids.join(",")),
                    Some(regulatory_features.iter().join(",")),
                    Some(repeats.iter().join(",")),
                    Some(*num_split_reads),
                    Some(*coverage),
                    breakpoint_sequence.clone(),
                ),
                SegmentAnnotation::Junction {
                    num_split_reads,
                    breakpoint_sequence,
                } => (
                    Some("split".to_string()),
                    None,
                    None,
                    None,
                    None,
                    None,
                    Some(*num_split_reads),
                    None,
                    Some(breakpoint_sequence.clone()),
                ),
            }
        } else {
            (None, None, None, None, None, None, None, None, None)
        };
        let is_coverage_segment = if let Some(k) = segment_type.as_ref() {
            k == "coverage"
        } else {
            false
        };
        let entry = SegmentEntry {
            graph_id: annotated_circle.graph_id(),
            circle_id: annotated_circle.circle_id(),
            kind: segment_type,
            target_from: reference.tid_to_tname(*from_ref_id).into(),
            from: *from,
            target_to: reference.tid_to_tname(*to_ref_id).into(),
            to: *to,
            length: if is_coverage_segment {
                Some(from.abs_diff(*to))
            } else {
                None
            },
            num_exons,
            gene_ids,
            gene_names,
            regulatory_features,
            repeats,
            coverage,
            num_split_reads,
            breakpoint_sequence,
        };
        writer.serialize(entry)?;
    }
    Ok(())
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
#[derive(Serialize, Debug, FieldNamesAsArray)]
struct FlatCircleTableInfo {
    event_id: String,
    graph_id: usize,
    circle_id: usize,
    circle_length: u32,
    segment_count: usize,
    regions: String,
    num_exons: usize,
    gene_ids: String,
    gene_names: String,
    regulatory_features: String,
    repeats: String,
    num_split_reads: usize,
    varlociraptor_event_prob: f32,
    prob_absent: f32,
    prob_artifact: f32,
    af_nanopore: Option<f32>,
}
