use std::collections::HashMap;
use std::iter::once;

use bam::{Header, Record};
use bio::data_structures::interval_tree::ArrayBackedIntervalTree;
use itertools::Itertools;
use num_traits::{One, Zero};
use rayon::prelude::*;

pub(crate) type ReferenceId = u32;
pub(crate) type Count = i32;
pub(crate) type ReadDepth = HashMap<ReferenceId, Vec<Count>>;
pub(crate) type SplitReadStorage = HashMap<String, Vec<Record>>;
pub(crate) type RecordTree<'a> = ArrayBackedIntervalTree<usize, &'a Record>;
pub(crate) type RecordTrees<'a> = HashMap<ReferenceId, RecordTree<'a>>;

pub(crate) fn read_depths<I: IntoIterator<Item = Record>>(
    record_iter: I,
    header: &bam::Header,
    exact: bool,
) -> anyhow::Result<ReadDepth> {
    let mut coverage: HashMap<_, _> = header
        .reference_names()
        .iter()
        .zip(header.reference_lengths())
        .map(|(name, length)| {
            (
                header.reference_id(name).unwrap(),
                vec![Count::zero(); *length as usize],
            )
        })
        .collect();
    if exact {
        record_iter.into_iter().for_each(|record| {
            let ref_id = record.ref_id();
            if ref_id >= 0 {
                let counts = coverage
                    .get_mut(&(ref_id as u32))
                    .expect("no such reference id");
                for (_seq_pos, ref_pos) in record.matching_pairs() {
                    let s = &mut counts[ref_pos as usize];
                    *s = s.saturating_add(Count::one());
                }
            }
        });
    } else {
        record_iter.into_iter().for_each(|record| {
            let start = record.start();
            let end = record.calculate_end().saturating_sub(1);
            let ref_id = record.ref_id();
            if ref_id >= 0 && start >= 0 && end >= 0 {
                assert!(start < end);
                let counts = coverage
                    .get_mut(&(ref_id as u32))
                    .expect("no such reference id");
                let s = &mut counts[start as usize];
                *s = s.saturating_add(Count::one());
                let e = &mut counts[end as usize];
                *e = e.saturating_sub(Count::one());
            }
        });
        coverage.values_mut().for_each(|counts| {
            counts.iter_mut().fold(Count::zero(), |acc, c| {
                *c = c.saturating_add(acc);
                *c
            });
        });
    }
    Ok(coverage)
}

pub(crate) struct SplitReadInfo<'a> {
    pub(crate) split_reads: &'a HashMap<String, Vec<Record>>,
    pub(crate) record_trees: RecordTrees<'a>,
    name_to_id: HashMap<String, u32>,
}

impl<'a> SplitReadInfo<'a> {
    pub(crate) fn new(header: &Header, split_reads: &'a SplitReadStorage) -> anyhow::Result<Self> {
        let name_to_id = (0..header.n_references() as u32)
            .map(|i| (header.reference_name(i).unwrap().to_string(), i))
            .collect();

        let mut s = Self {
            split_reads,
            record_trees: HashMap::new(),
            name_to_id,
        };

        // Iterate all records in the iterator and build both
        // an interval tree for start..end -> records lookups
        // as well as coverage information, both stratified by reference/chrom/target
        let trees = &mut s.record_trees;
        eprintln!("Building split read info");
        let records = s.split_reads.values().flatten();
        for record in records {
            let ref_id = record.ref_id();
            if ref_id < 0 {
                continue;
            }
            let (start, end) = (record.start(), record.calculate_end());
            if start >= 0 && end >= 0 {
                let (start, end) = (start as usize, end as usize);
                if end <= start {
                    dbg!(&record);
                    continue;
                }
                let tree = trees.entry(ref_id as u32).or_default();
                tree.insert(start..end, record);
            }
        }
        trees.values_mut().for_each(|tree| tree.index());

        Ok(Self {
            split_reads: s.split_reads,
            record_trees: s.record_trees,
            name_to_id: s.name_to_id,
        })
    }

    #[allow(unused)]
    pub(crate) fn record_tree(&self, contig: &str) -> &RecordTree {
        &self.record_trees[&self.name_to_id[contig]]
    }

    pub(crate) fn coverage_regions(
        read_depths: &ReadDepth,
        min_read_depth: usize,
    ) -> HashMap<u32, Vec<u32>> {
        read_depths
            .par_iter()
            .map(|(ref_id, counts)| {
                let groups = counts
                    .iter()
                    .enumerate()
                    .group_by(|&(_pos, depth)| (*depth as usize) < min_read_depth);
                (
                    *ref_id,
                    groups
                        .into_iter()
                        .map(|(_below_min_read_depth, mut group)| {
                            let first_pos = group.next().unwrap().0.min(counts.len() - 1);
                            first_pos as u32 // - if _below_min_read_depth { 0 } else { 1 }
                        })
                        .chain(once((counts.len() - 1) as u32))
                        .collect_vec(),
                )
            })
            .collect()
    }
}
