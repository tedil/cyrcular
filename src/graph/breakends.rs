use std::collections::HashMap;
use std::io::BufWriter;
use std::path::PathBuf;

use clap::Parser;
use serde::Serialize;

use crate::graph::cli::{EventId, GraphStorage};
use crate::graph::{plot, Cycle, GraphCaller};
use crate::util::vcf_header_from_bam;

#[derive(Parser)]
pub struct BreakendArgs {
    /// Input BAM file
    #[clap(parse(from_os_str))]
    input: PathBuf,

    /// Input reference FASTA
    #[clap(long, parse(from_os_str))]
    reference: PathBuf,

    /// Output BCF
    #[clap(long, parse(from_os_str))]
    output: PathBuf,

    /// Output graph in msgpack format
    #[clap(long, parse(from_os_str))]
    graph: PathBuf,

    /// Output connected components in dot format in this folder
    #[clap(long, parse(from_os_str))]
    dot: Option<PathBuf>,

    /// Number of reads used to distinguish between covered and not-covered regions.
    #[clap(long, default_value = "2")]
    min_read_depth: usize,

    /// Minimum number of split reads needed as evidence for potential breakpoint
    #[clap(long, default_value = "3")]
    min_split_reads: usize,

    /// Maximum number of plausible paths per component/subgraph
    #[clap(long, default_value = "5")]
    max_paths_per_component: usize,

    /// Maximum length of deletions between neighbouring breakpoints
    #[clap(long)]
    max_deletion_length: Option<usize>,

    #[clap(short, long, default_value = "0")]
    threads: u16,
}

pub fn main_breakends(args: BreakendArgs) -> anyhow::Result<()> {
    let records_reader = bam::IndexedReader::build()
        .additional_threads(args.threads)
        .from_path(&args.input)?;
    // let reference_reader = bio::io::fasta::IndexedReader::from_file(&args.reference)?;

    let vcf_path = args.output.clone();
    let graph_path = args.graph.clone();

    let mut graph_caller = GraphCaller {
        // reference: reference_reader,
        reference_path: args.reference.clone(),
        records: records_reader,
        min_read_depth: args.min_read_depth,
        min_split_reads: args.min_split_reads,
        max_paths_per_component: args.max_paths_per_component,
        max_deletion_length: args.max_deletion_length,
    };
    let graph = graph_caller.build_graph()?;

    eprintln!("Generating breakend events from graph");
    let events = graph_caller.breakends(&graph);
    eprintln!("Number of events: {}", events.len());

    let valid_paths: HashMap<EventId, Vec<Cycle>> = events
        .iter()
        .map(|event| (event.id, event.cycles.clone()))
        .collect();
    eprintln!(
        "Number of valid paths: {}",
        &valid_paths.values().map(|v| v.len()).sum::<usize>()
    );

    let mut path_serializer = std::fs::File::create(&graph_path)
        .map(BufWriter::new)
        .map(rmp_serde::Serializer::new)?;
    GraphStorage {
        // graph,
        valid_paths,
    }
    .serialize(&mut path_serializer)?;

    eprintln!("Writing vcf");
    let mut vcf_writer = std::fs::File::create(vcf_path)
        .map(BufWriter::new)
        .map(noodles::bcf::Writer::new)?;
    let bam_header = graph_caller.records.header().clone();
    let header = vcf_header_from_bam(&bam_header);
    vcf_writer.write_file_format()?;
    vcf_writer.write_header(&header)?;
    let raw_header = header.to_string();
    let string_map = raw_header.parse()?;
    events.iter().for_each(|event| {
        event
            .breakpoint_records
            .iter()
            .flatten()
            .for_each(|record| {
                vcf_writer
                    .write_vcf_record(&header, &string_map, record)
                    .expect("failed writing bcf record");
            })
    });

    if let Some(dot_dir) = &args.dot {
        eprintln!("Writing dot files");
        std::fs::create_dir_all(dot_dir)?;
        for (event_id, component) in events.into_iter().map(|event| (event.id, event.subgraph)) {
            let dot = plot::graph_to_dot(&component, &bam_header);
            let path = dot_dir.join(format!("{event_id}.dot", event_id = event_id));
            std::fs::write(path, dot)?;
        }
    }
    Ok(())
}
