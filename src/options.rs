//! Code for supporting program configuration/settings.
use clap::ArgMatches;

/// Options for the "coverage" command.
#[derive(Clone, Debug)]
pub struct Options {
    // Generic arguments
    /// Number of background threads to use for compression/decompression in I/O.
    pub io_threads: u32,

    // I/O related
    /// Path to sample BAM file.
    pub input: Vec<String>,
    /// Path to coverage depth BCF file.
    pub output: String,
    /// Path to PED file (for creating `FOUNDER_*` counts).
    pub input_ped: Option<String>,
    /// Path to panel information file (for creating by-population counts), format is
    /// "SAMPLE<tab>SUB-POPULATION<tab>POPULATION<ignored>"
    pub input_panel: Option<String>,
    /// Path to input FAi file.
    pub input_fasta: Option<String>,
}

impl Options {
    /// Build options from ArgMatches.
    pub fn new(matches: &ArgMatches) -> Self {
        Self {
            io_threads: matches
                .value_of("io_threads")
                .unwrap()
                .parse::<u32>()
                .unwrap(),

            input: matches
                .values_of("input")
                .expect("Problem getting input paths from command line")
                .map(|s| s.to_string())
                .collect::<Vec<_>>(),
            output: matches
                .value_of("output")
                .expect("Problem getting output path from command line")
                .to_string(),
            input_ped: matches.value_of("input_ped").map(|x| x.to_string()),
            input_panel: matches.value_of("input_panel").map(|x| x.to_string()),
            input_fasta: matches.value_of("input_fasta").map(|x| x.to_string()),
        }
    }
}
