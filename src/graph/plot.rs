use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

use anyhow::Result;
use bam::{Header, IndexedReader};
use clap::Parser;
use petgraph::dot::Dot;

use crate::graph::cli::{EventId, GraphStorage};
use crate::graph::{BreakpointGraph, Cycle, EdgeType};
use crate::plot::detailed_coverage;

#[derive(Parser)]
pub struct PlotArgs {
    /// Input BAM file
    #[clap(parse(from_os_str))]
    input: PathBuf,

    /// Input graph / valid paths
    #[clap(long, parse(from_os_str))]
    graph: PathBuf,

    /// Output directory
    #[clap(long, parse(from_os_str))]
    output: PathBuf,

    #[clap(short, long, default_value = "0")]
    threads: u16,
}

pub fn main_plot(args: PlotArgs) -> anyhow::Result<()> {
    let graph: GraphStorage =
        rmp_serde::from_read(std::fs::File::open(&args.graph).map(BufReader::new)?)?;
    let mut records_reader = bam::IndexedReader::build()
        .additional_threads(args.threads)
        .from_path(&args.input)?;
    let header = records_reader.header().clone();
    std::fs::create_dir_all(&args.output)?;

    for (event_id, cycles) in graph.valid_paths {
        for cycle in cycles {
            let divs = plot_breakend(&args, &mut records_reader, &header, event_id, &cycle)?;
            let header = r##"<html><head><meta charset="utf-8"/></head><body>
            <script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/MathJax.js?config=TeX-AMS-MML_SVG"></script>
            <script type="text/javascript">if (window.MathJax) {MathJax.Hub.Config({SVG: {font: "STIX-Web"}});}</script>
            <script type="text/javascript">window.PlotlyConfig = {MathJaxConfig: 'local'};</script>
            <script src="https://cdn.plot.ly/plotly-1.54.6.min.js"></script>"##;
            let divs = divs.join("\n");
            let html = format!("{header}{divs}</body></html>", header = header, divs = divs);
            let output = args.output.join(format!(
                "graph_{event_id}_{cycle_id}.html",
                event_id = event_id,
                cycle_id = cycle.id
            ));
            std::fs::write(output, html)?;
        }
    }
    Ok(())
}

fn plot_breakend(
    args: &PlotArgs,
    records_reader: &mut IndexedReader<File>,
    header: &Header,
    event_id: EventId,
    cycle: &Cycle,
) -> Result<Vec<String>> {
    let mut divs = Vec::new();
    for (i, &((from_id, from_pos), (to_id, to_pos), edge)) in cycle.edges.iter().enumerate() {
        if edge.edge_type.contains(EdgeType::Neighbour) {
            // skip edges across chromosome borders
            // as well as the last edge which is always a backwards edge
            if from_id != to_id || i == cycle.edges.len() - 1 {
                continue;
            }
            let (from_pos, to_pos) = if from_pos > to_pos {
                (to_pos, from_pos)
            } else {
                (from_pos, to_pos)
            };
            let from_name = header
                .reference_name(from_id)
                .ok_or_else(|| anyhow::anyhow!("no such ref_id"))?;
            let ref_len = header
                .reference_len(from_id)
                .ok_or_else(|| anyhow::anyhow!("no such ref_id"))?;
            let original_region = crate::plot::Region::new(from_name.into(), from_pos, to_pos);
            let plot_region = original_region.clone();
            // bam_region might be larger if flank > 0, todo
            let bam_region = bam::Region::new(from_id, from_pos, to_pos);
            let bin_size = auto_bin_size(&original_region);
            let mut detailed_coverage =
                detailed_coverage(records_reader, &plot_region, &bam_region, ref_len, bin_size)?;
            let plot = crate::plot::plot(
                &mut detailed_coverage,
                &original_region,
                &plot_region,
                bin_size,
                0,
            )?;
            let plot_id = format!(
                "graph_{event_id}_{cycle_id}_{edge}",
                event_id = event_id,
                cycle_id = cycle.id,
                edge = i
            );
            let filename = format!("{plot_id}.html", plot_id = plot_id,);
            let output = args.output.join(filename);
            let div = plot.to_inline_html(None);
            divs.push(div);
            plot.to_html(output);
        }
    }
    Ok(divs)
}

fn auto_bin_size(region: &crate::plot::Region) -> usize {
    let mut bin_size = (region.len() as f64).sqrt().round() as usize;
    let max_bins = 1200;
    if region.len() / bin_size > max_bins {
        bin_size = region.len() / max_bins;
    }
    bin_size
}

pub fn graph_to_dot(graph: &BreakpointGraph) -> String {
    let dot = Dot::with_attr_getters(
        &graph,
        &[],
        &|_graph, (_, _, edge_weight)| {
            let t = &edge_weight.edge_type;
            format!(
                "color=\"{}\"",
                if t.contains(EdgeType::Neighbour) && t.contains(EdgeType::Split) {
                    "#46A473"
                } else if t.contains(EdgeType::Neighbour) && t.contains(EdgeType::Deletion) {
                    "#CB7459"
                } else if !t.contains(EdgeType::Neighbour) && t.contains(EdgeType::Split) {
                    "#A38F2D"
                } else if t.contains(EdgeType::Neighbour)
                    && !t.contains(EdgeType::Split)
                    && !t.contains(EdgeType::Deletion)
                {
                    "#7E87D6"
                } else {
                    "black"
                }
            )
        },
        &|_, _| String::new(),
    );
    format!("{:?}", dot)
}
