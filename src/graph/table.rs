use std::collections::HashMap;
use std::fs::File;
use clap::Parser;
use serde::{Deserialize, Serialize};
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use anyhow::Result;
use itertools::Itertools;
use csv;
use crate::cli::{CircleId, EventId};
use crate::common::ReferenceId;
use crate::graph::annotate::{AnnotatedCircle, GraphId, SegmentAnnotation, VarlociraptorInfo};
use crate::graph::Position;

#[derive(Parser)]
struct TableArgs {
    /// Path for the overview circle table
    #[clap(long, parse(from_os_str))]
    circle_table: PathBuf,

    /// Output directory for detailed per-segment information for each circle
    #[clap(long, parse(from_os_str))]
    segment_tables: PathBuf,
}

fn tables_from_annotated_graph() {
    //
    // eprintln!("Reading breakend records");
    // let (mut reader, header, string_maps) = read_bcf(&args)?;
    // // FIXME: this assumes that the order of contigs is consistent between the vcf header, the graph, the bam and other files.
    // let tid_to_tname = tid_to_tname(&header);
    //
    // let event_records = group_event_records(&mut reader, &header, &string_maps);
    //
    // std::fs::create_dir_all(&args.segment_tables)?;
    // eprintln!("Building circle table");
    // let circle_table = graph
    //     .valid_paths
    //     .into_iter()
    //     .flat_map(|(graph_id, circles)| {
    //         circles
    //             .into_iter()
    //             .enumerate()
    //             .filter_map(|(circle_id, circle)| {
    //                 let event_name = format!("graph_{}_circle_{}", graph_id, circle_id);
    //                 if let Some(records) = event_records.get(&event_name) {
    //                     let entry = circle_table_entry(
    //                         &tid_to_tname,
    //                         &gene_annotations,
    //                         &regulatory_annotations,
    //                         &reference,
    //                         graph_id,
    //                         circle_id,
    //                         circle,
    //                         records,
    //                         args.breakpoint_sequence_length,
    //                     );
    //                     let mut segment_writer = segment_table_writer(&args, &event_name)
    //                         .unwrap_or_else(|_| panic!("Failed creating file for {}", &event_name));
    //                     write_segment_table_into(&mut segment_writer, &entry, &tid_to_tname)
    //                         .unwrap();
    //                     Some(entry)
    //                 } else {
    //                     eprintln!("WARNING: Event {} not found in breakend VCF", event_name);
    //                     None
    //                 }
    //             })
    //             .collect_vec()
    //     });
    //
    // let circle_table = circle_table
    //     .into_iter()
    //     .map(|(graph_id, circle_id, circle_info, segment_annotations)| {
    //         let (num_exons, gene_ids, gene_names, num_split_reads, regions, regulatory_features) =
    //             collapse_segment_annotations(&tid_to_tname, circle_info, &segment_annotations);
    //         let gene_ids = gene_ids.into_iter().unique().join(",");
    //         let gene_names = gene_names.into_iter().unique().join(",");
    //         let regulatory_features = regulatory_features.into_iter().join(",");
    //         let regions = if circle_info.segment_count == 1 {
    //             regions.get(0).cloned().unwrap_or_default()
    //         } else {
    //             regions.join(",")
    //         };
    //         FlatCircleTableInfo {
    //             event_id: format!("{}-{}", graph_id, circle_id),
    //             graph_id,
    //             circle_id,
    //             circle_length: circle_info.length,
    //             segment_count: circle_info.segment_count,
    //             regions,
    //             score: circle_info.score,
    //             num_split_reads,
    //             num_exons,
    //             gene_ids,
    //             gene_names,
    //             regulatory_features,
    //             prob_present: circle_info.prob_present,
    //             prob_absent: circle_info.prob_absent,
    //             prob_artifact: circle_info.prob_artifact,
    //             af_nanopore: circle_info.af_nanopore,
    //         }
    //     })
    //     .sorted_unstable_by_key(|table_entry| {
    //         (
    //             OrderedFloat(1. - table_entry.prob_present),
    //             -(table_entry.num_exons as i64),
    //             -(table_entry.score as i64),
    //             OrderedFloat(table_entry.prob_absent),
    //             table_entry.graph_id,
    //             table_entry.circle_id,
    //         )
    //     })
    //     .collect_vec();
    // write_circle_table(
    //     File::create(&args.circle_table).map(BufWriter::new)?,
    //     &circle_table,
    // )?;

    // let gene_ids = gene_ids.into_iter().unique().join(",");
    // let gene_names = gene_names.into_iter().unique().join(",");
    // let regulatory_features = regulatory_features.into_iter().join(",");
    // let regions = if annotated_circle.annotated_paths.len() <= 2 {
    //     regions.get(0).cloned().unwrap_or_default()
    // } else {
    //     regions.join(",")
    // };

    // let varlociraptor_annotations = varlociraptor_info(records);
}
fn segment_table_writer(
    args: &TableArgs,
    event_name: &String,
) -> Result<csv::Writer<BufWriter<File>>> {
    Ok(csv::WriterBuilder::new().delimiter(b'\t').from_writer(
        File::create(
            &args
                .segment_tables
                .join(format!("{}_segments.tsv", event_name)),
        )
        .map(BufWriter::new)?,
    ))
}

type CircleTableEntry = (
    EventId,
    VarlociraptorInfo,
    AnnotatedCircle,
);

fn write_circle_table<W: Write>(writer: W, circle_table: &[FlatCircleTableInfo]) -> Result<()> {
    let mut writer: csv::Writer<_> = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(writer);
    for entry in circle_table {
        writer.serialize(entry)?;
    }
    Ok(())
}

fn write_segment_table_into<W: Write>(
    writer: &mut csv::Writer<W>,
    circle_entry: &CircleTableEntry,
    tid_to_tname: &HashMap<ReferenceId, String>,
) -> Result<()> {
    let (_, circle_info, annotated_circle) = circle_entry;
    #[derive(Serialize, Debug)]
    struct SegmentEntry {
        graph_id: GraphId,
        circle_id: CircleId,
        circle_length: u32,
        segment_count: usize,
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
                    num_split_reads,
                    coverage,
                    breakpoint_sequence,
                } => (
                    Some("coverage".to_string()),
                    Some(*num_exons),
                    Some(gene_names.join(",")),
                    Some(gene_ids.join(",")),
                    Some(regulatory_features.iter().join(",")),
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
                    Some(*num_split_reads),
                    None,
                    Some(breakpoint_sequence.clone()),
                ),
            }
        } else {
            (None, None, None, None, None, None, None, None)
        };
        let is_coverage_segment = if let Some(k) = segment_type.as_ref() {
            k == "coverage"
        } else {
            false
        };
        let entry = SegmentEntry {
            graph_id: annotated_circle.graph_id(),
            circle_id: annotated_circle.circle_id(),
            circle_length: annotated_circle.circle_length(),
            segment_count: annotated_circle.num_segments(),
            kind: segment_type,
            target_from: tid_to_tname[from_ref_id].clone(),
            from: *from,
            target_to: tid_to_tname[to_ref_id].clone(),
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
#[derive(Serialize, Debug)]
struct FlatCircleTableInfo {
    event_id: String,
    graph_id: usize,
    circle_id: usize,
    circle_length: usize,
    segment_count: usize,
    regions: String,
    score: usize,
    num_exons: usize,
    gene_ids: String,
    gene_names: String,
    regulatory_features: String,
    num_split_reads: usize,
    prob_present: f32,
    prob_absent: f32,
    prob_artifact: f32,
    af_nanopore: Option<f32>,
}
