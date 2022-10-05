use std::collections::HashMap;
use std::fs::File;
use std::path::PathBuf;

use anyhow::Result;
use bam::{Header, IndexedReader};
use clap::Parser;
use itertools::Itertools;

use crate::graph::cli::{EventId, GraphStorage};
use crate::graph::{Breakpoint, BreakpointGraph, Cycle, Edge, EdgeInfo, EdgeType};
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

    #[clap(short, long, default_value = "3")]
    breakpoint_margin: u32,

    #[clap(short, long, default_value = "0")]
    threads: u16,
}

pub fn main_plot(args: PlotArgs) -> anyhow::Result<()> {
    let graph = GraphStorage::from_path(&args.graph)?;
    let mut records_reader = bam::IndexedReader::build()
        .additional_threads(args.threads)
        .from_path(&args.input)?;
    let header = records_reader.header().clone();
    std::fs::create_dir_all(&args.output)?;

    for (event_id, cycles) in graph.valid_paths {
        for cycle in cycles {
            let divs = plot_breakend(
                &args,
                &mut records_reader,
                &header,
                event_id,
                &cycle,
                args.breakpoint_margin,
            )?;
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
    breakpoint_margin: u32,
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
            let mut detailed_coverage = detailed_coverage(
                records_reader,
                &plot_region,
                &bam_region,
                ref_len,
                bin_size,
                breakpoint_margin,
            )?;
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

pub fn graph_to_dot(graph: &BreakpointGraph, header: &Header) -> String {
    let legend = r##"{
  rankdir=LR;
  rank=same;
  node [shape=plaintext]
  subgraph cluster_01 { 
    label = "Legend";
    key [label=<<table border="0" cellpadding="2" cellspacing="0" cellborder="0">
      <tr><td align="right" port="i1">neighbour</td></tr>
      <tr><td align="right" port="i2">split</td></tr>
      <tr><td align="right" port="i3">neighbour + split</td></tr>
      <tr><td align="right" port="i4">deletion</td></tr>
      </table>>]
    key2 [label=<<table border="0" cellpadding="2" cellspacing="0" cellborder="0">
      <tr><td port="i1">&nbsp;</td></tr>
      <tr><td port="i2">&nbsp;</td></tr>
      <tr><td port="i3">&nbsp;</td></tr>
      <tr><td port="i4">&nbsp;</td></tr>
      </table>>]
    key:i1:e -> key2:i1:w [color="#7E87D6"]
    key:i2:e -> key2:i2:w [color="#A38F2D", style=dotted]
    key:i3:e -> key2:i3:w [color="#46A473", style=dashed]
    key:i4:e -> key2:i4:w [color="#CB7459", style=dotted]
    }
}"##;

    // build node identifiers from reference id and positions for unique node names in dot repr.
    let node_id = |node: &(u32, u32)| -> String { format!("node_{}_{}", node.0, node.1) };

    // sort nodes by reference id and position, then generate node definitions
    let nodes = graph.nodes().sorted_unstable().map(|node| {
        let (ref_id, pos) = node;
        let label = format!(
            "{}:{}",
            header
                .reference_name(ref_id)
                .unwrap_or(&format!("[{}]", ref_id)),
            pos
        );
        format!("{} [label=\"{}\", shape = box]", node_id(&node), label)
    });

    // for edges, we also merge directed edges if both to -> from and from -> to edges exist
    // with the same properties.
    let mut edge_map: HashMap<(Breakpoint, Breakpoint), (Edge, bool)> = HashMap::new();

    // generate edge label from edge information
    let edge_label = |edge: &(Breakpoint, Breakpoint, EdgeInfo)| -> String {
        let (from, to, edge_info) = edge;
        let t = edge_info.edge_type;
        let (color, style) = if t.contains(EdgeType::Neighbour) && t.contains(EdgeType::Split) {
            ("#46A473", Some("dashed"))
        } else if t.contains(EdgeType::Neighbour) && t.contains(EdgeType::Deletion) {
            ("#CB7459", Some("dotted"))
        } else if !t.contains(EdgeType::Neighbour) && t.contains(EdgeType::Split) {
            ("#A38F2D", Some("dotted"))
        } else if t.contains(EdgeType::Neighbour)
            && !t.contains(EdgeType::Split)
            && !t.contains(EdgeType::Deletion)
        {
            ("#7E87D6", None)
        } else {
            ("black", None)
        };
        let label = format!("{:?}", edge_info);
        format!(
            "{} -> {} [label=\"{}\", color=\"{}\",{} fontsize = \"12\", minlen = \"{}\", penwidth = \"{:.2}\"]",
            node_id(from),
            node_id(to),
            label,
            color,
            style.map(|s| format!(" style=\"{}\",", s)).unwrap_or_default(),
            1.max(((edge_info.distance as f64 + 1.).log10() / 2.).round() as usize),
            1f64.max((edge_info.coverage + 1.).log(4.)),
        )
    };

    // if both from -> to and to -> from directed edges exist, merge into one edge
    graph
        .all_edges()
        .sorted_unstable_by_key(|(from, to, _)| (*from, *to))
        .for_each(|(from, to, edge_info)| {
            if let Some((rev_edge, bidirectional)) = edge_map.get_mut(&(to, from)) {
                if edge_info.distance == rev_edge.2.distance
                    && edge_info.num_split_reads == rev_edge.2.num_split_reads
                    && (edge_info.coverage - rev_edge.2.coverage).abs() < 1e-4
                {
                    *bidirectional = true;
                } else {
                    edge_map.insert((from, to), ((from, to, *edge_info), false));
                }
            } else {
                edge_map.insert((from, to), ((from, to, *edge_info), false));
            }
        });

    // generate edge dot representations, set arrowheads according to directionality
    let edges = edge_map.values().map(|(edge, bidirectional)| {
        if *bidirectional {
            format!("{} [dir=both]", edge_label(edge))
        } else {
            edge_label(edge)
        }
    });

    format!(
        r##"digraph {{
          graph [fontname = "helvetica"];
          node [fontname = "helvetica"];
          edge [fontname = "helvetica"];
          {legend}
          {{
            rank = same;
            rankdir = LR;
            {graph}
          }}
        }}"##,
        legend = legend,
        graph = nodes.chain(edges).join(";\n"),
    )
}
