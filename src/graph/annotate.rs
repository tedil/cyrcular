use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::ops::Index;
use std::path::{Path, PathBuf};

use anyhow::Result;
use bio::data_structures::interval_tree::ArrayBackedIntervalTree;
use bio::io::gff;
use bio::io::gff::GffType::GFF3;
use clap::Parser;
use flate2::bufread::MultiGzDecoder;
use itertools::Itertools;
use noodles::vcf::Record;
use serde::Deserialize;
use serde::Serialize;

use crate::cli::{CircleId, GraphStorage};
use crate::common::ReferenceId;
use crate::graph::EdgeType::{Neighbour, Split};
use crate::graph::{Cycle, Position};

#[derive(Parser)]
pub(crate) struct AnnotateArgs {
    /// Input graph in msgpack format
    #[clap(parse(from_os_str))]
    graph: PathBuf,

    /// Reference FASTA file
    #[clap(long, parse(from_os_str))]
    reference: PathBuf,

    /// (b)gzipped GFF3 file containing gene annotations with respect to the reference sequence
    #[clap(long, parse(from_os_str))]
    gene_annotation: PathBuf,

    /// (b)gzipped GFF3 file containing regulatory annotations with respect to the reference sequence
    #[clap(long, parse(from_os_str))]
    regulatory_annotation: PathBuf,

    /// Length of the sequence flanking the breakpoint
    #[clap(long, default_value = "2000")]
    breakpoint_sequence_length: u32,

    /// VCF/BCF file containing filtered and processed breakend events for the circles described in the graph
    #[clap(long, parse(from_os_str))]
    breakends: Option<PathBuf>,

    #[clap(long)]
    output: PathBuf,
}

pub(crate) fn main_annotate(args: AnnotateArgs) -> Result<()> {
    eprintln!("Reading graph");
    let graph = GraphStorage::from_path(&args.graph)?;

    eprintln!("Reading reference");
    let reference = Reference::from_path(&args.reference)?;

    eprintln!("Reading annotation");
    let gene_annotations = read_gff3(&args.gene_annotation, |record| {
        record.feature_type() == "gene" || record.feature_type() == "exon"
    })?;
    let regulatory_annotations = read_gff3(&args.regulatory_annotation, |_| true)?;
    let annotated_graph = annotate_graph(
        graph,
        &reference,
        &gene_annotations,
        &regulatory_annotations,
        None,
        args.breakpoint_sequence_length as usize,
    );
    annotated_graph.to_path(&args.output)?;
    Ok(())
}

fn annotate_graph(
    graph: GraphStorage,
    reference: &Reference,
    gene_annotations: &Annotations,
    regulatory_annotations: &Annotations,
    event_grouped_records: Option<HashMap<String, Vec<Record>>>,
    breakpoint_sequence_length: usize,
) -> AnnotatedGraph {
    if let Some(_event_records) = event_grouped_records {
        todo!();
    }

    eprintln!("Building circle table");
    let annotated_paths = graph
        .valid_paths
        .into_iter()
        .map(|(graph_id, circles)| {
            let annotated_circles = circles
                .into_iter()
                .map(|circle| {
                    let circle_id = circle.id;
                    let segment_annotations = segment_annotation(
                        gene_annotations,
                        regulatory_annotations,
                        reference,
                        circle,
                        breakpoint_sequence_length as u32,
                    );
                    AnnotatedCircle {
                        graph_id,
                        circle_id,
                        annotated_paths: segment_annotations,
                    }
                })
                .collect_vec();
            (graph_id, annotated_circles)
        })
        .collect_vec();
    let annotated_graph = AnnotatedGraph {
        annotated_circles: annotated_paths,
    };
    annotated_graph
}

pub(crate) type GraphId = usize;

#[derive(Debug, Serialize, Deserialize, Clone)]
pub(crate) struct AnnotatedCircle {
    graph_id: GraphId,
    circle_id: CircleId,
    annotated_paths: Vec<(Segment, Option<SegmentAnnotation>)>,
}

impl AnnotatedCircle {
    pub(crate) fn circle_id(&self) -> CircleId {
        self.circle_id
    }

    pub(crate) fn graph_id(&self) -> CircleId {
        self.graph_id
    }
}

impl AnnotatedCircle {
    pub(crate) fn summary(&self, reference: &Reference) -> CircleSummary {
        let num_segments = self.num_segments();
        let (num_exons, gene_ids, gene_names, num_split_reads, regions, regulatory_features) =
            self.annotated_paths.iter().fold(
                (0, vec![], vec![], 0, vec![], HashSet::new()),
                |(mut n_ex, mut gids, mut gnames, mut n_splits, mut regions, mut r_feats),
                 (segment, sa)| {
                    match sa {
                        Some(SegmentAnnotation::Segment {
                            num_exons,
                            gene_ids,
                            gene_names,
                            regulatory_features,
                            num_split_reads,
                            coverage: _,
                            breakpoint_sequence: _,
                        }) => {
                            n_ex += num_exons;
                            gids.extend(gene_ids.clone());
                            gnames.extend(gene_names.clone());
                            n_splits += if num_segments == 1 {
                                // if there's only "one" segment, the split reads loop from one end of the segment to its other end
                                *num_split_reads
                            } else {
                                // otherwise, split reads are described by junction edges and counted below
                                0
                            };
                            r_feats.extend(regulatory_features.clone());
                            let (tid_from, from, tid_to, to) = segment;
                            assert_eq!(tid_from, tid_to);
                            regions.push(format!(
                                "{}:{}-{}",
                                reference.tid_to_tname(*tid_from),
                                from,
                                to
                            ));
                        }
                        Some(SegmentAnnotation::Junction {
                            num_split_reads, ..
                        }) => {
                            n_splits += num_split_reads;
                        }
                        _ => (),
                    }
                    (n_ex, gids, gnames, n_splits, regions, r_feats)
                },
            );
        CircleSummary {
            event_id: event_id(self.graph_id, self.circle_id),
            graph_id: self.graph_id,
            circle_id: self.circle_id,
            circle_length: self.circle_length(),
            num_exons,
            gene_ids,
            gene_names,
            num_split_reads,
            regions,
            regulatory_features,
            segment_count: self.num_segments(),
        }
    }

    pub(crate) fn segments(&self) -> impl Iterator<Item = &(Segment, Option<SegmentAnnotation>)> {
        self.annotated_paths.iter()
    }
}

pub(crate) fn event_id(graph_id: GraphId, circle_id: CircleId) -> String {
    format!("{}-{}", graph_id, circle_id)
}

impl AnnotatedCircle {
    pub(crate) fn circle_length(&self) -> u32 {
        self.annotated_paths
            .iter()
            .map(|(segment, segment_annotation)| {
                if matches!(segment_annotation, Some(SegmentAnnotation::Segment { .. })) {
                    segment_length(segment).unwrap_or(0)
                } else {
                    0
                }
            })
            .sum()
    }

    pub(crate) fn num_segments(&self) -> usize {
        let n_segs = self
            .annotated_paths
            .iter()
            .filter(|(_, segment_annotation)| {
                matches!(segment_annotation, Some(SegmentAnnotation::Segment { .. }))
            })
            .count();
        if n_segs == self.annotated_paths.len() && n_segs == 2 {
            // if there are no junctions, then the circle is a single segment
            1
        } else {
            n_segs
        }
    }
}

#[derive(Debug, Serialize, Deserialize)]
pub(crate) struct AnnotatedGraph {
    pub(crate) annotated_circles: Vec<(GraphId, Vec<AnnotatedCircle>)>,
}

impl AnnotatedGraph {
    pub(crate) fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        Ok(rmp_serde::from_read(
            std::fs::File::open(path).map(BufReader::new)?,
        )?)
    }

    pub(crate) fn to_path<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let mut serializer = std::fs::File::create(&path)
            .map(BufWriter::new)
            .map(rmp_serde::Serializer::new)?;
        self.serialize(&mut serializer)?;
        Ok(())
    }
}

pub(crate) struct CircleSummary {
    pub(crate) event_id: String,
    pub(crate) graph_id: GraphId,
    pub(crate) circle_id: CircleId,
    pub(crate) circle_length: u32,
    pub(crate) segment_count: usize,
    pub(crate) regions: Vec<String>,
    pub(crate) num_split_reads: NumSplitReads,
    pub(crate) num_exons: NumExons,
    pub(crate) gene_ids: Vec<String>,
    pub(crate) gene_names: Vec<String>,
    pub(crate) regulatory_features: HashSet<String>,
}

impl AnnotatedGraph {
    pub(crate) fn summaries(&self, reference: &Reference) -> Vec<CircleSummary> {
        self.annotated_circles
            .iter()
            .flat_map(|(_graph_id, segment_annotations)| {
                segment_annotations
                    .iter()
                    .map(|annotated_circle| annotated_circle.summary(reference))
                    .collect_vec()
            })
            .collect_vec()
    }
}

type Annotation = ArrayBackedIntervalTree<u64, gff::Record>;
type Annotations = HashMap<String, Annotation>;

fn read_gff3<P: AsRef<Path> + std::fmt::Debug>(
    path: P,
    filter: fn(&gff::Record) -> bool,
) -> Result<Annotations> {
    let mut annotations = gff::Reader::new(
        File::open(path)
            .map(BufReader::new)
            .map(MultiGzDecoder::new)
            .unwrap(),
        GFF3,
    );
    let mut trees: HashMap<String, Annotation> = HashMap::new();
    annotations
        .records()
        .map(|r| r.unwrap())
        .filter(filter)
        .for_each(|r| {
            let tree = trees
                .entry(r.seqname().to_string())
                .or_insert_with(Annotation::new);
            let (start, end) = (r.start().min(r.end()), r.start().max(r.end()));
            tree.insert(*start..*end, r);
        });
    trees.values_mut().for_each(|tree| tree.index());
    Ok(trees)
}

pub(crate) struct Reference {
    inner: HashMap<ReferenceId, Vec<u8>>,
    tid_to_tname: HashMap<ReferenceId, String>,
    tname_to_tid: HashMap<String, ReferenceId>,
}

impl Reference {
    pub(crate) fn tid_to_tname(&self, tid: ReferenceId) -> &str {
        &self.tid_to_tname[&tid]
    }
}

impl Index<ReferenceId> for Reference {
    type Output = [u8];
    fn index(&self, index: ReferenceId) -> &Self::Output {
        &self.inner[&index]
    }
}

impl Index<&str> for Reference {
    type Output = [u8];
    fn index(&self, index: &str) -> &Self::Output {
        &self.inner[&self.tname_to_tid[index]]
    }
}

impl Reference {
    pub(crate) fn from_path<P: AsRef<Path> + std::fmt::Debug>(path: P) -> Result<Self> {
        let inner: HashMap<ReferenceId, (String, Vec<u8>)> =
            bio::io::fasta::Reader::from_file(path)?
                .records()
                .enumerate()
                .map(|(i, r)| {
                    r.map(|r| (i as u32, (r.id().to_string(), r.seq().to_vec())))
                        .map_err(|e| e.into())
                })
                .collect::<Result<_>>()?;
        let tid_to_tname: HashMap<ReferenceId, String> = inner
            .iter()
            .map(|(tid, (tname, _))| (*tid, tname.clone()))
            .collect();
        let tname_to_tid: HashMap<String, ReferenceId> = tid_to_tname
            .iter()
            .map(|(tid, tname)| (tname.clone(), *tid))
            .collect();
        let inner = inner
            .into_iter()
            .map(|(tid, (_, seq))| (tid, seq))
            .collect();
        Ok(Self {
            inner,
            tid_to_tname,
            tname_to_tid,
        })
    }
}

type NumExons = usize;
type NumSplitReads = usize;

type Segment = (ReferenceId, Position, ReferenceId, Position);

fn segment_length(segment: &Segment) -> Option<u32> {
    let (tid_from, from, tid_to, to) = segment;
    (tid_from == tid_to).then(|| to.abs_diff(*from))
}

fn segment_annotation(
    gene_annotations: &Annotations,
    regulatory_annotations: &Annotations,
    reference: &Reference,
    circle: Cycle,
    breakpoint_sequence_length: u32,
) -> Vec<(Segment, Option<SegmentAnnotation>)> {
    circle
        .edges
        .iter()
        .map(|&((from_ref_id, from), (to_ref_id, to), edge_info)| {
            let is_coverage_segment = edge_info.edge_type.contains(Neighbour)
                && edge_info.coverage > 1e-4
                && from_ref_id == to_ref_id;
            if is_coverage_segment {
                let breakpoint_sequence = (to < from).then(|| {
                    breakpoint_sequence(
                        reference,
                        breakpoint_sequence_length,
                        from_ref_id,
                        from,
                        to_ref_id,
                        to,
                    )
                    .unwrap()
                });

                let chrom = reference.tid_to_tname(from_ref_id);
                let gene_annotation = gene_annotations
                    .get(chrom)
                    .or_else(|| gene_annotations.get(&format!("chr{}", chrom)));
                let regulatory_annotation = regulatory_annotations
                    .get(chrom)
                    .or_else(|| regulatory_annotations.get(&format!("chr{}", chrom)));
                let (exons, gene_ids, gene_names) = if let Some(gene_annotation) = gene_annotation {
                    annotate_genes_and_exons(from, to, gene_annotation)
                } else {
                    (HashSet::default(), vec![], vec![])
                };
                let regulatory_features = if let Some(regulatory_annotation) = regulatory_annotation
                {
                    annotate_regulatory_features(from, to, regulatory_annotation)
                } else {
                    HashSet::default()
                };

                (
                    (from_ref_id, from, to_ref_id, to),
                    Some(SegmentAnnotation::Segment {
                        num_exons: exons.len(),
                        gene_ids,
                        gene_names,
                        regulatory_features,
                        coverage: edge_info.coverage,
                        num_split_reads: edge_info.num_split_reads,
                        breakpoint_sequence,
                    }),
                )
            } else {
                // not a coverage segment but a split junction
                if edge_info.edge_type.contains(Split) {
                    let breakpoint_sequence = breakpoint_sequence(
                        reference,
                        breakpoint_sequence_length,
                        from_ref_id,
                        from,
                        to_ref_id,
                        to,
                    )
                    .unwrap();
                    let annot = SegmentAnnotation::Junction {
                        num_split_reads: edge_info.num_split_reads,
                        breakpoint_sequence,
                    };
                    ((from_ref_id, from, to_ref_id, to), Some(annot))
                } else {
                    ((from_ref_id, from, to_ref_id, to), None)
                }
            }
        })
        .sorted_unstable_by_key(|((from_ref_id, from, to_ref_id, to), _)| {
            (*from_ref_id, *from, *to_ref_id, *to)
        })
        .collect_vec()
}

fn annotate_genes_and_exons(
    from: Position,
    to: Position,
    annot: &Annotation,
) -> (HashSet<String>, Vec<String>, Vec<String>) {
    let (from2, to2) = (from.min(to), from.max(to));
    let (exons, gene_ids, gene_names) = annot
        .find(from2 as u64..to2 as u64)
        .iter()
        .map(|entry| entry.data())
        .filter_map(|data| match data.feature_type() {
            "gene" => Some((
                None,
                data.attributes().get("gene_id").cloned(),
                data.attributes().get("Name").cloned(),
            )),
            "exon" => Some((data.attributes().get("exon_id").cloned(), None, None)),
            _ => None,
        })
        .fold(
            (HashSet::new(), vec![], vec![]),
            |(mut exons, mut gene_ids, mut gene_names), (exon_id, gene_id, gene_name)| {
                if let Some(exon_id) = exon_id {
                    exons.insert(exon_id);
                };
                gene_ids.extend(gene_id);
                gene_names.extend(gene_name);
                (exons, gene_ids, gene_names)
            },
        );
    let gene_ids = gene_ids.into_iter().unique().collect_vec();
    let gene_names = gene_names.into_iter().unique().collect_vec();
    (exons, gene_ids, gene_names)
}

fn annotate_regulatory_features(
    from: Position,
    to: Position,
    annot: &Annotation,
) -> HashSet<String> {
    let (from2, to2) = (from.min(to), from.max(to));
    annot
        .find(from2 as u64..to2 as u64)
        .iter()
        .map(|entry| entry.data().feature_type().to_string())
        .collect()
}

fn breakpoint_sequence(
    reference: &Reference,
    breakpoint_sequence_length: u32,
    from_ref_id: ReferenceId,
    from: Position,
    to_ref_id: ReferenceId,
    to: Position,
) -> Result<String> {
    let breakpoint_sequence_length = breakpoint_sequence_length / 2;
    let from_seq = &reference[from_ref_id];
    let seq1 = &from_seq[from.saturating_sub(breakpoint_sequence_length) as usize..from as usize];
    let to_seq = &reference[to_ref_id];
    let to_len = to_seq.len();
    let seq2 = &to_seq
        [to as usize..((to.saturating_add(breakpoint_sequence_length)) as usize).min(to_len)];
    Ok(String::from_utf8([seq1, seq2].concat())?)
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub(crate) enum SegmentAnnotation {
    Segment {
        num_exons: usize,
        gene_ids: Vec<String>,
        gene_names: Vec<String>,
        regulatory_features: HashSet<String>,
        num_split_reads: usize,
        coverage: f64,
        breakpoint_sequence: Option<String>,
    },
    Junction {
        num_split_reads: usize,
        breakpoint_sequence: String,
    },
}
