use std::cmp::Ordering;
use std::collections::HashMap;
use std::fmt::{Debug, Formatter};
use std::hash::Hash;
use std::io::{Read, Seek};
use std::iter::from_fn;
use std::iter::FromIterator;
use std::path::PathBuf;
use std::time::Instant;

use anyhow::Result;
use bam::Record;
use bio::data_structures::interval_tree::ArrayBackedIntervalTree;
use enumflags2::{bitflags, BitFlags};
use indexmap::set::IndexSet;
use itertools::Itertools;
use lazy_static::lazy_static;
use noodles::vcf::header::info::Key;
use noodles::vcf::record::info::field::Value::Integer;
use ordered_float::OrderedFloat;
use petgraph::prelude::*;
use petgraph::visit::{Data, IntoNeighborsDirected, NodeCount};
use priority_queue::double_priority_queue::DoublePriorityQueue;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::common::{read_depths, ReadDepth, ReferenceId, SplitReadInfo};
use crate::util::{default_filter, split_reads};

mod annotate;
mod breakends;
pub(crate) mod cli;
mod plot;
mod table;

type Position = u32;
type Breakpoint = (ReferenceId, Position);

struct GraphCaller<R: Seek + Read> {
    // reference: bio::io::fasta::IndexedReader<R>,
    reference_path: PathBuf,
    records: bam::IndexedReader<R>,
    min_read_depth: usize,
    min_split_reads: usize,
    max_paths_per_component: usize,
    max_deletion_length: Option<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Cycle {
    id: usize,
    edges: Vec<Edge>,
}

#[derive(Debug)]
pub struct BreakendEvent {
    id: usize,
    subgraph: BreakpointGraph,
    cycles: Vec<Cycle>,
    breakpoint_records: Vec<Vec<noodles::vcf::Record>>,
}

type Support = usize;

impl<R: Seek + Read> GraphCaller<R> {
    pub(crate) fn build_graph(&mut self) -> Result<BreakpointGraph> {
        let records_reader = &mut self.records;
        let header = records_reader.header().clone();
        let now = Instant::now();
        let read_depths = read_depths(
            records_reader.full_by(default_filter).flatten(),
            &header,
            false,
        )?;
        let n_pos = read_depths.values().map(|c| c.len() as f64).sum::<f64>();
        let mean_read_depth = read_depths
            .values()
            .flatten()
            .map(|c| *c as f64)
            .sum::<f64>()
            / n_pos;

        eprintln!("mean read depth: {:.4}", mean_read_depth);
        eprintln!("Calculating coverage took: {:#?}", now.elapsed());

        let now = Instant::now();
        let split_reads = split_reads(records_reader)?;
        eprintln!("Found {} split reads", &split_reads.len());
        // TODO: move split_reads into SplitReadInfo
        let split_read_info = SplitReadInfo::new(records_reader.header(), &split_reads)?;

        eprintln!("Building split-read info took: {:#?}", now.elapsed());

        let graph = build_graph(
            &split_read_info,
            &read_depths,
            self.min_read_depth,
            self.min_split_reads,
            self.max_deletion_length,
        );
        Ok(graph)
    }

    pub(crate) fn breakends(&mut self, graph: &BreakpointGraph) -> Vec<BreakendEvent> {
        let is_valid_path = |path: &Vec<Breakpoint>| -> bool {
            use EdgeType::*;

            let edge_weights = path
                .windows(2)
                .map(|w| *graph.edge_weight(w[0], w[1]).unwrap())
                .collect_vec();

            let mut last_edge: Option<BitFlags<EdgeType>> = Option::None;
            let valid = edge_weights.iter().all(|edge| {
                let this_edge = edge.edge_type;
                let valid = if let Some(last_edge) = last_edge {
                    last_edge.contains(Neighbour)
                        && (this_edge.contains(Split) | this_edge.contains(Deletion))
                        || ((last_edge.contains(Split) || last_edge.contains(Deletion))
                            && this_edge.contains(Neighbour))
                } else {
                    true
                };

                last_edge = Some(this_edge);
                valid
            });

            let has_some_coverage = edge_weights.iter().step_by(2).all(|weight| {
                weight.edge_type.contains(Neighbour) && weight.distance > 0 && weight.coverage > 0.
            }) || edge_weights.iter().skip(1).step_by(2).all(|weight| {
                weight.edge_type.contains(Neighbour) && weight.distance > 0 && weight.coverage > 0.
            });
            valid && has_some_coverage
        };

        let mut components = petgraph::algo::kosaraju_scc(&graph);
        eprintln!("Number of components: {}", components.len());
        let reference_path = self.reference_path.clone();
        let header = self.records.header().clone();
        let max_paths = self.max_paths_per_component;

        components
            .par_iter_mut()
            .enumerate()
            .filter(|(_, component)| component.len() < 100)
            .map(|(event_id, component)| {
                component.sort_unstable();
                (event_id, component)
            })
            .map(|(event_id, component)| (event_id, crate::graph::subgraph(graph, component)))
            .filter(|(_, subgraph)| subgraph.edge_count() < 100 && subgraph.node_count() < 100)
            .map(|(event_id, subgraph): (usize, BreakpointGraph)| {
                // Generate all cyclic paths starting and ending at the same breakpoint node.
                // While generating paths, order children such that edges with more support are
                // walked earlier than those with less support.
                let start_nodes = subgraph.nodes().filter(|node| {
                    let num_incoming = subgraph.edges_directed(*node, Direction::Incoming).count();
                    let num_outgoing = subgraph.edges_directed(*node, Direction::Outgoing).count();
                    !(num_incoming == 1 && num_outgoing == 1) || subgraph.node_count() == 2
                });
                let valid_paths = start_nodes
                    .flat_map(|node| {
                        all_cycles_from(&subgraph, node, |graph, parent, child_a, child_b| {
                            // while walking paths, order children before pushing them onto the stack
                            // such that children more likely to yield valid & 'good' cycles are visited
                            // first
                            let edge_a = graph.edge_weight(parent, *child_a).unwrap();
                            let edge_b = graph.edge_weight(parent, *child_b).unwrap();
                            OrderedFloat(edge_a.coverage)
                                .cmp(&OrderedFloat(edge_b.coverage))
                                .then(edge_a.num_split_reads.cmp(&edge_b.num_split_reads))
                                .reverse()
                        })
                    })
                    .filter(is_valid_path);

                type BreakpointPath = Vec<Breakpoint>;

                let mut pq = DoublePriorityQueue::<BreakpointPath, Support>::new();
                let score = |path: &BreakpointPath| -> Support {
                    let s = path
                        .windows(2)
                        .map(|w| (w[0], w[1], *graph.edge_weight(w[0], w[1]).unwrap()))
                        .map(|(_a, _b, edge)| {
                            let read_depth_support = (edge.coverage + 1.).log10();
                            let split_read_support = (edge.num_split_reads as f64).sqrt();
                            read_depth_support + split_read_support + 1.
                        })
                        .product::<f64>()
                        .powf(1. / path.len() as f64);
                    s.round() as usize
                };

                for path in valid_paths {
                    let s = score(&path);
                    // while the priority queue isn't full, simply insert
                    if pq.len() < max_paths {
                        pq.push(path, s);
                        continue;
                    }
                    // check if the current one is better than the worst
                    // if so, ditch the worst one for the current one
                    // otherwise, stop.
                    if let Some((_, worst_score)) = pq.peek_min() {
                        if s > *worst_score {
                            pq.pop_min();
                        } else {
                            break;
                        }
                    }
                    pq.push(path, s);
                }

                let (records, cycles): (Vec<_>, Vec<_>) = pq
                    .into_iter()
                    .enumerate()
                    .map(|(path_id, (path, score))| {
                        (
                            path_id,
                            score,
                            path.windows(2)
                                .map(|w| (w[0], w[1], *graph.edge_weight(w[0], w[1]).unwrap()))
                                .collect_vec(),
                        )
                    })
                    .filter(|(_, _, edges)| {
                        edges.iter().any(|(_, _, w)| {
                            w.edge_type.contains(EdgeType::Neighbour)
                                && w.distance > 0
                                && w.coverage > 0.
                        })
                    })
                    .filter_map(|(path_id, score, path)| {
                        let mut reference =
                            bio::io::fasta::IndexedReader::from_file(&reference_path)
                                .expect("Failed opening fasta reader");
                        let mut records = breakend_event(
                            event_id,
                            path_id,
                            &path,
                            score,
                            &mut reference,
                            &header,
                        )
                        .unwrap_or_else(|_| {
                            panic!("Failed building breakend event for {}", event_id)
                        });
                        // TODO: find out why these circle with no proper coverage-segments
                        // do not get filtered out earlier
                        let has_length = records.iter().any(|record| {
                            if let Some(Integer(length)) = (*record.info())
                                .get(&*CIRCLE_LENGTH_KEY)
                                .and_then(|v| v.value())
                            {
                                *length > 0
                            } else {
                                false
                            }
                        });
                        if !has_length {
                            None
                        } else {
                            records.sort_unstable_by_key(|record| {
                                (
                                    record.chromosome().to_string(),
                                    usize::from(record.position()),
                                )
                            });
                            let cycle = Cycle {
                                id: path_id,
                                edges: path,
                            };
                            Some((records, cycle))
                        }
                    })
                    .unzip();
                BreakendEvent {
                    id: event_id,
                    subgraph,
                    cycles,
                    breakpoint_records: records,
                }
            })
            .collect()
    }
}

fn subgraph(graph: &BreakpointGraph, component: &[Breakpoint]) -> BreakpointGraph {
    let mut subgraph = graph.clone();
    let mut remove_nodes = Vec::new();
    for node in subgraph.nodes() {
        if !component.contains(&node) {
            remove_nodes.push(node);
        }
    }
    for node in remove_nodes {
        subgraph.remove_node(node);
    }
    subgraph
}

#[repr(u8)]
#[bitflags]
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
enum EdgeType {
    None = 0b0001,
    Neighbour = 0b0010,
    Deletion = 0b0100,
    Split = 0b1000,
}

#[derive(Copy, Clone, PartialEq, Serialize, Deserialize)]
pub struct EdgeInfo {
    coverage: f64,
    num_split_reads: usize,
    edge_type: BitFlags<EdgeType>,
    distance: u32,
}

impl Debug for EdgeInfo {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        f.write_fmt(format_args!(
            "cov: {:.2}, splits: {}, dist: {}",
            self.coverage, self.num_split_reads, self.distance
        ))
    }
}

type Edge = (Breakpoint, Breakpoint, EdgeInfo);

lazy_static! {
    pub(crate) static ref SUPPORT_KEY: Key = Key::Other("Support".into(),);
    pub(crate) static ref CIRCLE_LENGTH_KEY: Key = Key::Other("CircleLength".into(),);
    pub(crate) static ref CIRCLE_SEGMENT_COUNT_KEY: Key = Key::Other("CircleSegmentCount".into(),);
    pub(crate) static ref NUM_SPLIT_READS_KEY: Key = Key::Other("SplitReads".into(),);
}

fn breakend_event<R: Read + Seek>(
    event_id: usize,
    path_id: usize,
    path: &[Edge],
    score: Support,
    reference: &mut bio::io::fasta::IndexedReader<R>,
    bam_header: &bam::Header,
) -> Result<Vec<noodles::vcf::Record>> {
    use noodles::vcf;
    use noodles::vcf::record::info::field::*;
    use noodles::vcf::record::info::Field;
    use noodles::vcf::record::*;

    let mut ref_base_at = |ref_name: &str, idx: u64| -> char {
        reference
            .fetch(ref_name, idx, idx + 1)
            .unwrap_or_else(|_| panic!("Failed fetching {}:{}", ref_name, idx));
        let mut text = vec![0; 1];
        reference
            .read(&mut text)
            .unwrap_or_else(|_| panic!("Failed reading {}:{}", ref_name, idx));
        (text[0] as char).to_ascii_uppercase()
    };

    // common information for all sub-events on the path
    let event_id = format!("graph_{}_circle_{}", event_id, path_id);

    let mut events = Vec::new();

    // for each edge in the path that's not just a neighbour edge,
    // build two breakend events which are mates.
    // odd/even edges are covered segments, even/odd edges are splits/deletions,
    // so only keep every other edge.
    let edges = path;
    let first_is_neighbour = edges[0].2.edge_type.contains(EdgeType::Neighbour);

    // length of circle = sum of lengths of segments with coverage
    let (num_segments, circle_length) = edges
        .iter()
        .enumerate()
        .filter(|(i, _)| i % 2 != (first_is_neighbour as usize))
        .map(|(_, (_start, _stop, weight))| {
            if weight.edge_type.contains(EdgeType::Neighbour) {
                (1, weight.distance)
            } else {
                (0, 0)
            }
        })
        .fold((0u32, 0u32), |(n_segs, circ_len), (n, l)| {
            (n_segs + n, circ_len + l)
        });

    let edges = edges
        .iter()
        .enumerate()
        .filter(|(i, _)| i % 2 == (first_is_neighbour as usize))
        .map(|(_, e)| e)
        .collect_vec();

    let num_edges = edges.len();
    for (i, (from, to, edge)) in edges
        .into_iter()
        .filter(|(_, _, edge)| {
            num_edges == 1
                || edge.edge_type.contains(EdgeType::Split)
                || edge.edge_type.contains(EdgeType::Deletion)
        })
        .enumerate()
    {
        let partners = [*from, *to].map(|(ref_id, pos)| {
            let ref_name = bam_header.reference_name(ref_id).unwrap();
            let ref_base = ref_base_at(ref_name, pos as u64);
            let chrom_from = Chromosome::Name(ref_name.into());
            (chrom_from, ref_name, ref_base, pos)
        });
        let num_split_reads = if edge.edge_type.contains(EdgeType::Split) {
            edge.num_split_reads
        } else {
            0
        };
        let sub_event_id = i * partners.len();
        for (j, (chrom, _ref_name, ref_base, ref_pos)) in partners.iter().enumerate() {
            let infos = vec![
                Field::new((*SUPPORT_KEY).clone(), Some(Integer(score as i32))),
                Field::new(
                    (*CIRCLE_LENGTH_KEY).clone(),
                    Some(Integer(circle_length as i32)),
                ),
                Field::new(
                    (*CIRCLE_SEGMENT_COUNT_KEY).clone(),
                    Some(Integer(num_segments as i32)),
                ),
                Field::new(
                    (*NUM_SPLIT_READS_KEY).clone(),
                    Some(Integer(num_split_reads as i32)),
                ),
                Field::new(Key::SvType, Some(Value::String("BND".into()))),
                Field::new(Key::BreakendEventId, Some(Value::String(event_id.clone()))),
                Field::new(
                    Key::MateBreakendIds,
                    Some(Value::String(format!(
                        "{}_{}",
                        event_id,
                        sub_event_id + (1 - j)
                    ))),
                ),
            ];
            let builder = vcf::Record::builder()
                .set_quality_score(QualityScore::try_from(60.)?)
                .set_filters(Filters::Pass);

            let builder = {
                builder
                    .set_chromosome(chrom.clone())
                    .set_ids(format!("{}_{}", event_id, sub_event_id + j).parse()?)
                    // TODO we have to add 1 here since we're building **VCF** records which start indexing at 1, not 0
                    .set_position(Position::try_from(*ref_pos as usize + 1)?)
                    .set_reference_bases(format!("{}", ref_base).parse()?)
                    .set_info(Info::try_from(infos)?)
            };
            let [(_chrom_from, ref_name_from, _ref_base_from, ref_pos_from), (_chrom_to, ref_name_to, _ref_base_to, ref_pos_to)] =
                &partners;
            let builder = if j == 1 {
                // TODO we have to add 1 here since we're building **VCF** records which start indexing at 1, not 0
                builder.set_alternate_bases(
                    format!("]{}:{}]{}", ref_name_from, ref_pos_from + 1, ref_base).parse()?,
                )
            } else {
                // TODO we have to add 1 here since we're building **VCF** records which start indexing at 1, not 0
                builder.set_alternate_bases(
                    format!("{}[{}:{}[", ref_base, ref_name_to, ref_pos_to + 1).parse()?,
                )
            };
            events.push(builder.build()?);
        }
    }
    Ok(events)
}

type BreakpointGraph = DiGraphMap<Breakpoint, EdgeInfo>;

fn build_graph(
    split_info: &SplitReadInfo,
    read_depth: &ReadDepth,
    min_read_depth: usize,
    min_split_reads: usize,
    max_deletion_length: Option<usize>,
) -> BreakpointGraph {
    // Step 1: Determine breakpoints by coverage
    eprintln!("Step 1: Determining breakpoints (by coverage)");
    let mut now = Instant::now();
    let breakpoints = SplitReadInfo::coverage_regions(read_depth, min_read_depth)
        .into_iter()
        .flat_map(|(ref_id, transients)| transients.into_iter().map(move |pos| (ref_id, pos)))
        .unique()
        .collect_vec();

    let mut graph = BreakpointGraph::new();
    let mut breakpoints = breakpoints.to_vec();
    breakpoints.par_sort_unstable();
    eprintln!("{:#?}", now.elapsed());

    // Step 2: For each breakpoint, insert a corresponding node into the graph
    eprintln!("Step 2: Insert breakpoint nodes");
    eprintln!("(number of breakpoints: {})", breakpoints.len());
    now = Instant::now();
    for breakpoint in &breakpoints {
        graph.add_node(*breakpoint);
    }
    eprintln!("{:#?}", now.elapsed());

    // Step 3: For each pair of neighbouring breakpoints, insert an Edge
    eprintln!("Step 3: Inserting edges between neighbouring breakpoints");
    now = Instant::now();
    insert_neighbour_edges(
        min_read_depth,
        &mut graph,
        &breakpoints,
        read_depth,
        max_deletion_length,
    );
    eprintln!("{:#?}", now.elapsed());

    // Step 4: For each breakpoint, check if there are split reads which either start or end there;
    // and, if so, add edges between the corresponding nodes (adding those, if necessary)
    eprintln!("Step 4: Inserting edges between split-read nodes");
    now = Instant::now();
    insert_split_edges(split_info, read_depth, &mut graph, &mut breakpoints);
    eprintln!("{:#?}", now.elapsed());

    eprintln!("Step 5: Pruning graph, first removing certain edges, then removing nodes with no edges left (or not part of a circle)");
    now = Instant::now();
    let graph = prune(graph, min_split_reads);
    eprintln!("{:#?}", now.elapsed());
    graph
}

fn insert_split_edges(
    split_info: &SplitReadInfo,
    read_depth: &ReadDepth,
    graph: &mut BreakpointGraph,
    breakpoints: &mut [Breakpoint],
) {
    let mut empty_split_tree = ArrayBackedIntervalTree::new();
    empty_split_tree.index();
    let mut empty_bp_tree = ArrayBackedIntervalTree::new();
    empty_bp_tree.index();
    let margin = 15;
    let breakpoints_by_ref = breakpoints
        .iter()
        .copied()
        .into_group_map_by(|(ref_id, _)| *ref_id);
    let breakpoint_trees: HashMap<_, _> = breakpoints_by_ref
        .into_iter()
        .map(|(ref_id, breakpoints)| {
            (
                ref_id,
                breakpoints
                    .iter()
                    .map(|(_, pos)| ((*pos..*pos + 1), ()))
                    .collect::<ArrayBackedIntervalTree<u32, ()>>(),
            )
        })
        .collect();

    breakpoints.iter().for_each(|breakpoint| {
        let (ref_id, pos) = *breakpoint;
        let breakpoint_tree = breakpoint_trees.get(&ref_id).unwrap_or(&empty_bp_tree);
        let split_tree = split_info
            .record_trees
            .get(&ref_id)
            .unwrap_or(&empty_split_tree);
        let splits = split_tree
            .find(pos.saturating_sub(margin) as usize..pos.saturating_add(margin) as usize)
            .iter()
            .map(|entry| *entry.data())
            .collect_vec();
        let read_partners = read_partners(split_info, splits);
        for partner_read in &read_partners {
            let other_ref_id = partner_read.ref_id();
            let other_start_pos = partner_read.start();
            let other_end_pos = partner_read.calculate_end();
            assert!(other_ref_id >= 0 && other_start_pos >= 0 && other_end_pos >= 0);
            let possible_partners = [
                (other_ref_id as u32, other_start_pos as u32),
                (
                    other_ref_id as u32,
                    (other_end_pos as u32).saturating_sub(1),
                ),
            ];
            let mut possible_partners = possible_partners
                .iter()
                .filter(|partner| *partner != breakpoint)
                .filter_map(|partner| {
                    if !graph.contains_node(*partner) {
                        // try finding a node which is "close enough" first:
                        let (p_rid, p_pos) = *partner;
                        let close_nodes = breakpoint_tree
                            .find(p_pos.saturating_sub(margin)..p_pos.saturating_add(margin));
                        if !close_nodes.is_empty() {
                            // if there are multiple closest (unlikely), choose the one in the middle
                            Some((p_rid, close_nodes[close_nodes.len() / 2].interval().start))
                        } else {
                            // if there isn't, skip
                            None
                        }
                    } else {
                        Some(*partner)
                    }
                })
                .collect_vec();
            possible_partners.sort_unstable();
            for partner in possible_partners {
                if let Some(edge_info) = graph.edge_weight_mut(*breakpoint, partner) {
                    edge_info.edge_type |= EdgeType::Split;
                    edge_info.num_split_reads += 1;
                } else {
                    let coverage_mean = if partner.0 == ref_id {
                        neighbour_coverage(read_depth, ref_id, pos, partner.1)
                    } else {
                        0.
                    };
                    let mut edge_type = BitFlags::from_flag(EdgeType::Split);
                    if let Some(edge_info) = graph.edge_weight(partner, *breakpoint) {
                        if edge_info.edge_type.contains(EdgeType::Neighbour) {
                            edge_type |= EdgeType::Neighbour;
                        }
                    }
                    graph.add_edge(
                        *breakpoint,
                        partner,
                        EdgeInfo {
                            edge_type,
                            coverage: coverage_mean,
                            num_split_reads: 1,
                            distance: if ref_id == other_ref_id as u32 {
                                partner.1.abs_diff(pos) as u32
                            } else {
                                0
                            },
                        },
                    );
                }
            }
        }
    });
}

fn neighbour_coverage(read_depth: &ReadDepth, target: ReferenceId, a: u32, b: u32) -> f64 {
    let (a, b) = (a.min(b), a.max(b));
    let between_coverage = &read_depth[&target][a as usize..b as usize];
    let n = between_coverage.len() as f64;
    between_coverage
        .iter()
        .map(|count| (*count as f64) / n)
        .sum()
}

fn insert_neighbour_edges(
    min_read_depth: usize,
    graph: &mut GraphMap<Breakpoint, EdgeInfo, Directed>,
    breakpoints: &[(u32, u32)],
    read_depth: &ReadDepth,
    max_deletion_length: Option<usize>,
) {
    for window in breakpoints.windows(2) {
        let (&ba @ (target_a, pos_a), &bb @ (target_b, pos_b)) = (&window[0], &window[1]);
        // eprintln!("{:?} → {:?} ({})", &ba, &bb, bb.1 - ba.1);
        if target_a != target_b {
            continue;
        }
        let distance = pos_b - pos_a;
        let coverage_mean = neighbour_coverage(read_depth, target_a, pos_a, pos_b);
        let edge_type: BitFlags<EdgeType> = if coverage_mean < min_read_depth as f64 {
            if let Some(max_dist) = max_deletion_length {
                if distance as usize > max_dist {
                    continue;
                }
            }
            EdgeType::Neighbour | EdgeType::Deletion
        } else {
            BitFlags::from_flag(EdgeType::Neighbour)
        };
        graph.add_edge(
            ba,
            bb,
            EdgeInfo {
                coverage: coverage_mean,
                edge_type,
                num_split_reads: 0,
                distance,
            },
        );
    }
}

fn prune(
    mut graph: GraphMap<Breakpoint, EdgeInfo, Directed>,
    min_split_reads: usize,
) -> GraphMap<Breakpoint, EdgeInfo, Directed> {
    let components = petgraph::algo::kosaraju_scc(&graph);
    eprintln!("Number of components (pre-pruning): {}", components.len());

    // remove all edges with too few split-reads
    let to_remove = graph
        .all_edges()
        .filter_map(|(from, to, edge_info)| {
            if edge_info.num_split_reads > 0
                && edge_info.num_split_reads < min_split_reads
                && !edge_info.edge_type.contains(EdgeType::Deletion)
            {
                Some((from, to))
            } else {
                None
            }
        })
        .collect_vec();
    for (a, b) in &to_remove {
        graph.remove_edge(*a, *b);
    }
    let components = petgraph::algo::kosaraju_scc(&graph);
    eprintln!(
        "Number of components (after edge removal): {}",
        components.len()
    );

    // Remove all nodes that have no edges or are not part of a cycle
    let to_remove = graph
        .nodes()
        .par_bridge()
        .filter(|node| {
            let outgoing = graph.edges_directed(*node, Outgoing).collect_vec();
            if outgoing
                .iter()
                .any(|(_, _, edge)| edge.edge_type == EdgeType::Neighbour)
            {
                false
            } else {
                let incoming = graph.edges_directed(*node, Incoming).count();
                let too_few_edges = (incoming + outgoing.len()) <= 1;
                // let in_cycle = is_node_in_cycle(&graph, *node);
                let split_edges_only = split_edge_filter(&graph, node);
                too_few_edges || split_edges_only // || !in_cycle
            }
        })
        .collect::<Vec<_>>();
    for node in &to_remove {
        graph.remove_node(*node);
    }

    let components = petgraph::algo::kosaraju_scc(&graph);
    eprintln!(
        "Number of components (after node removal #1): {}",
        components.len()
    );

    let components = petgraph::algo::kosaraju_scc(&graph);
    let singletons = std::sync::atomic::AtomicUsize::new(0);
    let remove_nodes = components
        .par_iter()
        .flat_map(|component| {
            if component.len() <= 1 {
                singletons.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                component.clone()
            } else {
                component
                    .iter()
                    .filter(|&node| {
                        // remove nodes with split edges only
                        split_edge_filter(&graph, node)
                    })
                    .copied()
                    .collect_vec()
            }
        })
        .collect::<Vec<_>>();
    for node in remove_nodes {
        graph.remove_node(node);
    }

    let components = petgraph::algo::kosaraju_scc(&graph);
    eprintln!(
        "Number of components (after node removal #2): {}",
        components.len()
    );

    // Remove loops (self edges)
    let to_remove = graph
        .nodes()
        .filter_map(|node| {
            if graph.contains_edge(node, node) {
                Some((node, node))
            } else {
                None
            }
        })
        .collect_vec();
    for edge in to_remove {
        graph.remove_edge(edge.0, edge.1);
    }

    graph
}

fn split_edge_filter(graph: &BreakpointGraph, node: &Breakpoint) -> bool {
    [Direction::Incoming, Direction::Outgoing]
        .iter()
        .all(|dir| {
            graph.edges_directed(*node, *dir).all(|(_, _, edge)| {
                edge.edge_type.contains(EdgeType::Split)
                    && !edge.edge_type.contains(EdgeType::Neighbour)
                    && !edge.edge_type.contains(EdgeType::Deletion)
            })
        })
}

fn read_partners<'a>(split_info: &'a SplitReadInfo, splits: Vec<&Record>) -> Vec<&'a Record> {
    let read_partners = splits
        .iter()
        .flat_map(|read| {
            let name = std::str::from_utf8(read.name()).unwrap().to_owned();
            // skip the read itself. Since record does not implement Eq,
            // we check for alignments being different
            split_info.split_reads[&name]
                .iter()
                .filter(move |other| !read.aligned_pairs().eq(other.aligned_pairs()))
        })
        .collect_vec();
    read_partners
}

pub fn all_cycles_from<TargetColl, G, C>(
    graph: G,
    from: G::NodeId,
    cmp: C,
) -> impl Iterator<Item = TargetColl>
where
    G: NodeCount,
    G: Data,
    G: IntoNeighborsDirected,
    G::NodeId: Eq + Hash + Copy + Debug,
    TargetColl: FromIterator<G::NodeId>,
    C: Fn(G, G::NodeId, &G::NodeId, &G::NodeId) -> Ordering,
{
    let max_length = 30;
    let min_length = 2;

    // list of visited nodes
    let mut visited: IndexSet<G::NodeId> = IndexSet::from_iter(Some(from));
    // list of childs of currently exploring path nodes,
    // last elem is list of childs of last visited node
    let mut stack = vec![graph
        .neighbors_directed(from, Outgoing)
        .sorted_by(|child_a, child_b| cmp(graph, from, child_a, child_b))];

    from_fn(move || {
        while let Some(children) = stack.last_mut() {
            if let Some(child) = children.next() {
                if visited.len() < max_length {
                    if child == from {
                        if visited.len() >= min_length {
                            let path = visited
                                .iter()
                                .cloned()
                                .chain(Some(from))
                                .collect::<TargetColl>();
                            return Some(path);
                        }
                    } else if !visited.contains(&child) {
                        visited.insert(child);
                        stack.push(
                            graph
                                .neighbors_directed(child, Outgoing)
                                .sorted_by(|child_a, child_b| cmp(graph, child, child_a, child_b)),
                        );
                    }
                } else {
                    if (child == from || children.any(|v| v == from)) && visited.len() >= min_length
                    {
                        let path = visited
                            .iter()
                            .cloned()
                            .chain(Some(from))
                            .collect::<TargetColl>();
                        return Some(path);
                    }
                    stack.pop();
                    visited.pop();
                }
            } else {
                stack.pop();
                visited.pop();
            }
        }
        None
    })
}
