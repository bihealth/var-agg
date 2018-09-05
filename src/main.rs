// `error_chain!` can recurse deeply.
#![recursion_limit = "1024"]

// TODO: get contigs from FAI file

use std::collections::HashMap;
use std::fs::File;
use std::result;
use std::sync::atomic::Ordering;
use std::sync::{atomic, Arc};

extern crate bio;
use bio::io::fasta;

extern crate csv;

#[macro_use]
extern crate error_chain;

#[macro_use]
extern crate clap;
use clap::{App, ArgMatches};

#[macro_use]
extern crate slog;
extern crate slog_async;
extern crate slog_term;

use slog::Drain;

extern crate rust_htslib;
use rust_htslib::bcf::{self, Read};

mod bcf_utils;
mod options;
use options::Options;
mod pedigree;
use pedigree::Pedigree;

mod errors {
    // Create the Error, ErrorKind, ResultExt, and Result types
    error_chain!{}
}

pub use errors::*;

/// Custom `slog` Drain logic
struct RuntimeLevelFilter<D> {
    drain: D,
    log_level: Arc<atomic::AtomicIsize>,
}

impl<D> Drain for RuntimeLevelFilter<D>
where
    D: Drain,
{
    type Ok = Option<D::Ok>;
    type Err = Option<D::Err>;

    fn log(
        &self,
        record: &slog::Record,
        values: &slog::OwnedKVList,
    ) -> result::Result<Self::Ok, Self::Err> {
        let current_level = match self.log_level.load(Ordering::Relaxed) {
            0 => slog::Level::Warning,
            1 => slog::Level::Info,
            _ => slog::Level::Trace,
        };

        if record.level().is_at_least(current_level) {
            self.drain.log(record, values).map(Some).map_err(Some)
        } else {
            Ok(None)
        }
    }
}

/// Build bcf::Writer with appropriate header.
fn build_bcf_writer(
    path: &String,
    pedigree: &Option<&Pedigree>,
    populations: &Option<&Vec<String>>,
    contigs: &Option<&Vec<(String, i32)>>,
    reader: &bcf::Reader,
) -> Result<bcf::Writer> {
    let mut header = bcf::header::Header::from_template_subset(&reader.header(), &[])
        .chain_err(|| "Could not create header with subset of samples")?;

    let mut lines = vec![
        // Information over all individuals
        "##INFO=<ID=AC_Hemi,Number=A,Type=Integer,Description=\"Hemizygous Count\">",
        "##INFO=<ID=AC_Het,Number=A,Type=Integer,Description=\"Heterozygous Count\">",
        "##INFO=<ID=AC_Hom,Number=A,Type=Integer,Description=\"Homozygous Count\">",
    ];

    if pedigree.is_some() {
        lines.append(&mut vec![
            // Information over founders only
            "##INFO=<ID=FOUNDER_AC,Number=A,Type=Integer,Description=\"Alternate Allele Count \
             (in founders only)\">",
            "##INFO=<ID=FOUNDER_AN,Number=A,Type=Integer,Description=\"Total Allele Count \
             (in founders only)\">",
            "##INFO=<ID=FOUNDER_Hemi,Number=A,Type=Integer,Description=\"Hemizygous Count (in \
             founders only)\">",
            "##INFO=<ID=FOUNDER_Het,Number=A,Type=Integer,Description=\"Heterozygous Count (in \
             founders only)\">",
            "##INFO=<ID=FOUNDER_Hom,Number=A,Type=Integer,Description=\"Homozygous Count (in \
             founders only)\">",
        ]);
    }

    let mut more = Vec::new();
    if let Some(pops) = populations {
        let mut pops_copy: Vec<String> = (*pops).clone();
        pops_copy.insert(0, "POPMAX".to_string());
        for ref pop in pops_copy {
            more.push(format!(
                "##INFO=<ID={}_AC,Number=A,Type=Integer,Description=\"Alternate Allele Count \
                 (in founders only)\">",
                pop
            ));
            more.push(format!(
                "##INFO=<ID={}_AN,Number=A,Type=Integer,Description=\"Total Allele Count \
                 (in founders only)\">",
                pop
            ));
            more.push(format!(
                "##INFO=<ID={}_Hemi,Number=A,Type=Integer,Description=\"Hemizygous Count (in \
                 founders only)\">",
                pop
            ));
            more.push(format!(
                "##INFO=<ID={}_Het,Number=A,Type=Integer,Description=\"Heterozygous Count (in \
                 founders only)\">",
                pop
            ));
            more.push(format!(
                "##INFO=<ID={}_Hom,Number=A,Type=Integer,Description=\"Homozygous Count (in \
                 founders only)\">",
                pop
            ));
        }
    }

    for line in lines {
        header.push_record(line.as_bytes());
    }
    for line in more {
        header.push_record(line.as_bytes());
    }

    // Replace contig lines
    if let Some(contigs) = contigs {
        for (ref name, length) in contigs.iter() {
            header.remove_contig(name.as_bytes());
            header.push_record(format!("##contig=<ID={},length={}>", name, length).as_bytes());
        }
    }

    let uncompressed = !path.ends_with(".bcf") && !path.ends_with(".vcf.gz");
    let vcf = path.ends_with(".vcf") || path.ends_with(".vcf.gz");

    bcf::Writer::from_path(&path, &header, uncompressed, vcf)
        .chain_err(|| "Could not open BCF file for writing")
}

fn run(matches: ArgMatches) -> Result<()> {
    // Atomic variable controlling logging level
    let log_level = Arc::new(atomic::AtomicIsize::new(1));

    // Perform slog setup
    let decorator = slog_term::TermDecorator::new().build();
    let drain = slog_term::FullFormat::new(decorator).build();
    let drain = RuntimeLevelFilter {
        drain: drain,
        log_level: log_level.clone(),
    }.fuse();
    let drain = slog_async::Async::new(drain).build().fuse();

    let logger = slog::Logger::root(drain, o!());

    // Switch log level
    if matches.is_present("quiet") {
        log_level.store(0, Ordering::Relaxed);
    } else {
        log_level.store(
            1 + matches.occurrences_of("verbose") as isize,
            Ordering::Relaxed,
        );
    };

    let options = Options::new(&matches);
    info!(logger, "Options: {:?}", &options);

    let pedigree = options.input_ped.as_ref().map_or(
        Ok(None),
        |path: &String| -> Result<Option<Pedigree>> {
            info!(logger, "Reading pedigree");
            let pedigree = Pedigree::from_path(&path)?;
            info!(
                logger,
                "Pedigree has {} members",
                pedigree.individuals.len()
            );
            Ok(Some(pedigree))
        },
    )?;

    let pop_map = options.input_panel.as_ref().map_or(
        Ok(None),
        |path: &String| -> Result<Option<HashMap<String, String>>> {
            info!(logger, "Reading panel map");
            let file =
                File::open(path).chain_err(|| format!("Could not open panel file: {}", &path))?;
            let mut reader = csv::ReaderBuilder::new()
                .delimiter(b'\t')
                .flexible(true)
                .has_headers(false)
                .from_reader(file);
            let mut pop_map = HashMap::new();
            for record in reader.records() {
                let record = record.chain_err(|| "Problem reading panel record from file")?;
                pop_map.insert(
                    record
                        .get(0)
                        .expect("Could not field 0 of record")
                        .to_string(),
                    record
                        .get(2)
                        .expect("Could not expect field 2 of record")
                        .to_string(),
                );
            }
            info!(logger, "Panel map has, {} entries", pop_map.len());
            Ok(Some(pop_map))
        },
    )?;
    let populations = pop_map.as_ref().map(|pop_map| {
        let mut populations = pop_map.values().map(|x| x.clone()).collect::<Vec<_>>();
        populations.sort();
        populations
    });

    let contigs: Option<Vec<(String, i32)>> = options.input_fasta.as_ref().map_or(
        Ok(None),
        |path: &String| -> Result<Option<Vec<(String, i32)>>> {
            info!(logger, "Reading FAI file");
            let ref_reader = fasta::IndexedReader::from_file(&path)
                .chain_err(|| "Loading FAI for genome regions failed")?;

            Ok(Some(
                ref_reader
                    .index
                    .sequences()
                    .iter()
                    .map(|ref seq| (seq.name.clone(), seq.len as i32))
                    .collect::<Vec<_>>(),
            ))
        },
    )?;

    info!(logger, "Opening output BCF/VCF file");
    let mut writer = {
        let reader = bcf::Reader::from_path(&options.input[0]).chain_err(|| {
            format!(
                "Could not open input file {} for reading",
                &options.input[0]
            )
        })?;
        build_bcf_writer(
            &options.output,
            &pedigree.as_ref(),
            &populations.as_ref(),
            &contigs.as_ref(),
            &reader,
        )?;
    };

    info!(logger, "Starting processing");

    bcf_utils::build_index(&logger, &options.output)?;
    info!(logger, "All done. Have a nice day!");

    Ok(())
}

fn main() {
    let yaml = load_yaml!("cli.yaml");
    let matches = App::from_yaml(yaml).get_matches();

    if let Err(ref e) = run(matches) {
        eprintln!("error: {}", e);

        for e in e.iter().skip(1) {
            eprintln!("caused by: {}", e);
        }

        // The backtrace is not always generated. Try to run this example
        // with `RUST_BACKTRACE=1`.
        if let Some(backtrace) = e.backtrace() {
            eprintln!("backtrace: {:?}", backtrace);
        }

        ::std::process::exit(1);
    }
}
