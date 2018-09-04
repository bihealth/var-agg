/// Code for supporting program configuration/settings.

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

            input: matches.values_of("input").expect("Problem getting input paths from command line").map(|s| s.to_string()).collect::<Vec<_>>(),
            output: matches.value_of("output").expect("Problem getting output path from command line").to_string(),
        }
    }
}
