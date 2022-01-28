use serde::{Deserialize, Serialize};
use std::collections::HashMap;

use crate::graph::breakends::BreakendArgs;
use crate::graph::plot::PlotArgs;
use crate::graph::{breakends, plot, Cycle};
use clap::Parser;

#[derive(Parser)]
pub(crate) struct GraphArgs {
    #[clap(subcommand)]
    subcommand: Command,
}

#[derive(Parser)]
enum Command {
    Breakends(BreakendArgs),
    Plot(PlotArgs),
}

pub(crate) fn main(args: GraphArgs) -> anyhow::Result<()> {
    match args.subcommand {
        Command::Breakends(args) => breakends::main_breakends(args),
        Command::Plot(args) => plot::main_plot(args),
    }
}

pub(crate) type EventId = usize;

#[derive(Serialize, Deserialize)]
pub(crate) struct GraphStorage {
    // graph: BreakpointGraph, // petgraph serde feature not complete, yet
    pub(crate) valid_paths: HashMap<EventId, Vec<Cycle>>,
}
