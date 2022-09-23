use std::collections::HashMap;
use std::io::BufReader;
use std::path::Path;

use anyhow::Result;
use clap::Parser;
use serde::{Deserialize, Serialize};

use crate::graph::annotate::AnnotateArgs;
use crate::graph::breakends::BreakendArgs;
use crate::graph::plot::PlotArgs;
use crate::graph::{annotate, breakends, plot, Cycle};

#[derive(Parser)]
pub(crate) struct GraphArgs {
    #[clap(subcommand)]
    subcommand: Command,
}

#[derive(Parser)]
enum Command {
    Breakends(BreakendArgs),
    Plot(PlotArgs),
    Annotate(AnnotateArgs),
}

pub(crate) fn main(args: GraphArgs) -> Result<()> {
    match args.subcommand {
        Command::Breakends(args) => breakends::main_breakends(args),
        Command::Plot(args) => plot::main_plot(args),
        Command::Annotate(args) => annotate::main_annotate(args),
    }
}

pub(crate) type EventId = usize;

#[derive(Serialize, Deserialize)]
pub(crate) struct GraphStorage {
    // graph: BreakpointGraph, // petgraph serde feature not complete, yet
    pub(crate) valid_paths: HashMap<EventId, Vec<Cycle>>,
}

impl GraphStorage {
    pub(crate) fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        Ok(rmp_serde::from_read(
            std::fs::File::open(path).map(BufReader::new)?,
        )?)
    }

    pub fn event(&self, id: EventId) -> Option<&Vec<Cycle>> {
        self.valid_paths.get(&id)
    }
}
