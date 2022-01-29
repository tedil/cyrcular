use anyhow::Result;
use clap::Parser;

use graph::cli::GraphArgs;
use plot::PlotArgs;

use crate::graph::cli;

mod common;
mod graph;
mod plot;
mod util;

#[derive(Parser)]
#[clap(about, version, author)]
struct Args {
    #[clap(subcommand)]
    subcommand: Subcommand,
}

#[derive(Parser)]
enum Subcommand {
    Plot(PlotArgs),
    Graph(GraphArgs),
}

fn main() -> Result<()> {
    let args: Args = Args::parse();
    match args.subcommand {
        Subcommand::Plot(p) => plot::main(p),
        Subcommand::Graph(g) => cli::main(g),
    }
}
