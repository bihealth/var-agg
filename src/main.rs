// `error_chain!` can recurse deeply.
#![recursion_limit = "1024"]

// TODO: get contigs from FAI file

use std::collections::HashMap;
use std::collections::HashSet;
use std::env;
use std::fs::File;
use std::mem;
use std::result;
use std::sync::atomic::Ordering;
use std::sync::{atomic, Arc};

extern crate bio;
use bio::io::fasta;

#[macro_use]
extern crate clap;
use clap::{App, ArgMatches};

extern crate csv;

#[macro_use]
extern crate error_chain;

extern crate ordered_float;
use ordered_float::OrderedFloat;

extern crate rust_htslib;
use rust_htslib::bcf::record::{Genotype, GenotypeAllele};
use rust_htslib::bcf::{self, Read};

#[macro_use]
extern crate slog;
extern crate slog_async;
extern crate slog_term;

use slog::Drain;

extern crate shlex;

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
    populations: &Option<&Vec<String>>,
    contigs: &Option<&Vec<(String, i32)>>,
    reader: &bcf::Reader,
) -> Result<bcf::Writer> {
    let mut header = bcf::header::Header::from_template_subset(&reader.header(), &[])
        .chain_err(|| "Could not create header with subset of samples")?;

    header.push_record(
        format!(
            "##varaggCommand={}",
            env::args()
                .map(|s| shlex::quote(&s).to_string())
                .collect::<Vec<String>>()
                .join(" ")
        ).as_bytes(),
    );

    let lines = vec![
        // The following line is only used to fix the bad 1000 genomes VCF files.
        "##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Genotype Likelihoods\">",
        // Information over founders (if PED is given) or all individuals (if not).
        "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Alternate Allele Count\">",
        "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total Allele Count\">",
        "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Alternate Allele Frequency\">",
        "##INFO=<ID=Hemi,Number=A,Type=Integer,Description=\"Hemizygous Count\">",
        "##INFO=<ID=Het,Number=A,Type=Integer,Description=\"Heterozygous Count\">",
        "##INFO=<ID=Hom,Number=A,Type=Integer,Description=\"Homozygous Count\">",
        "##INFO=<ID=POPMAX,Number=1,Type=String,Description=\"Population with highest alternative
         frequency (for multi-allelic sites, the sum of all allele frequencies is used)\">",
    ];

    let mut more = Vec::new();
    if let Some(pops) = populations {
        let mut pops_copy: Vec<String> = (*pops).clone();
        pops_copy.insert(0, "POPMAX".to_string());
        pops_copy.push("Oth".to_string());
        for ref pop in pops_copy {
            more.push(format!(
                "##INFO=<ID={}_AC,Number=A,Type=Integer,Description=\"Alternate Allele Count \
                 (in population {})\">",
                pop, pop
            ));
            more.push(format!(
                "##INFO=<ID={}_AN,Number=1,Type=Integer,Description=\"Total Allele Count \
                 (in population {})\">",
                pop, pop
            ));
            more.push(format!(
                "##INFO=<ID={}_AF,Number=A,Type=Float,Description=\"Alternate allele frequency \
                 (in population {})\">",
                pop, pop
            ));
            more.push(format!(
                "##INFO=<ID={}_Hemi,Number=A,Type=Integer,Description=\"Hemizygous Count (in
                 population {})\">",
                pop, pop
            ));
            more.push(format!(
                "##INFO=<ID={}_Het,Number=A,Type=Integer,Description=\"Heterozygous Count (in \
                 population {})\">",
                pop, pop
            ));
            more.push(format!(
                "##INFO=<ID={}_Hom,Number=A,Type=Integer,Description=\"Homozygous Count (in \
                 population {})\">",
                pop, pop
            ));
        }
    }

    let mut out_header = bcf::header::Header::new();
    for line in lines {
        out_header.push_record(line.as_bytes());
    }
    for line in more {
        out_header.push_record(line.as_bytes());
    }

    // Construct header with proper contigs.
    if let Some(contigs) = contigs {
        for (ref name, length) in contigs.iter() {
            out_header.push_record(format!("##contig=<ID={},length={}>", name, length).as_bytes());
        }
    }

    // Use trick for overriding contig and incorrect counts in output.
    unsafe {
        out_header.inner = rust_htslib::htslib::bcf_hdr_merge(out_header.inner, header.inner)
    };
    mem::swap(&mut out_header.inner, &mut header.inner);

    let uncompressed = !path.ends_with(".bcf") && !path.ends_with(".vcf.gz");
    let vcf = path.ends_with(".vcf") || path.ends_with(".vcf.gz");

    bcf::Writer::from_path(&path, &header, uncompressed, vcf)
        .chain_err(|| "Could not open BCF file for writing")
}

// TODO: this is insufficient as we need to do this per-alternative allele
/// Statistics for one group (could be population, founder, all).
#[derive(Debug)]
struct GroupStats {
    /// Name of the group
    pub name: String,
    /// Total observed heterozygous alleles
    pub het: Vec<i32>,
    /// Total observed homozygous alleles
    pub hom: Vec<i32>,
    /// Total observed hemizygous alleles
    pub hemi: Vec<i32>,
    /// Total allele count
    pub an: i32,
}

impl GroupStats {
    /// Create a new group only by name.
    ///
    /// ## Args
    ///
    /// - `usize` -- number of alternative alleles
    pub fn new(name: &str, length: usize) -> Self {
        return GroupStats {
            name: name.to_string(),
            het: vec![0; length],
            hom: vec![0; length],
            hemi: vec![0; length],
            an: 0,
        };
    }

    /// Return observed alternate allele count.
    pub fn ac(&self) -> Vec<i32> {
        (0..(self.het.len()))
            .map(|i| self.het[i] + 2 * self.hom[i] + self.hemi[i])
            .collect::<Vec<_>>()
    }

    /// Return allele frequencies.
    pub fn af(&self) -> Vec<f32> {
        self.ac()
            .iter()
            .map(|x| ((*x as f64) / (self.an as f64)) as f32)
            .collect::<Vec<_>>()
    }

    /// Register a `Genotype`.
    pub fn tally(&mut self, gt: &Genotype) -> () {
        assert!(gt.len() >= 1 && gt.len() <= 2);
        if gt.len() == 1 {
            // One allele only, could only be hemizygous.
            match gt[0] {
                GenotypeAllele::Unphased(i) | GenotypeAllele::Phased(i) => {
                    self.an += 1;
                    let i = i as usize;
                    if i > 0 {
                        self.hemi[i - 1] += 1;
                    }
                }
                _ => (),
            }
        } else {
            // gt.len() == 2
            // One allele only, could only be hemizygous.
            match (gt[0], gt[1]) {
                (GenotypeAllele::Unphased(i), GenotypeAllele::Unphased(j))
                | (GenotypeAllele::Phased(i), GenotypeAllele::Unphased(j))
                | (GenotypeAllele::Unphased(i), GenotypeAllele::Phased(j))
                | (GenotypeAllele::Phased(i), GenotypeAllele::Phased(j)) => {
                    // Case: both present
                    self.an += 2;
                    let i = i as usize;
                    let j = j as usize;
                    if i == 0 && j != 0 {
                        self.het[j - 1] += 1;
                    } else if i != 0 && j == 0 {
                        self.het[i - 1] += 1;
                    } else if i != 0 && j != 0 {
                        if i == j {
                            self.hom[i - 1] += 1;
                        } else {
                            self.het[i - 1] += 1;
                            self.het[j - 1] += 1;
                        }
                    }
                }
                (_, GenotypeAllele::Unphased(i))
                | (_, GenotypeAllele::Phased(i))
                | (GenotypeAllele::Unphased(i), _)
                | (GenotypeAllele::Phased(i), _) => {
                    // Case: one missing, one present
                    self.an += 1;
                    let i = i as usize;
                    if i != 0 {
                        self.het[i - 1] += 1;
                    }
                }
                _ => {
                    // Case: both missing -- nop
                }
            }
        }
    }
}

fn process_input(
    logger: &slog::Logger,
    path: &String,
    pedigree: &Option<&Pedigree>,
    pop_map: &Option<&HashMap<String, String>>,
    populations: &Option<&Vec<String>>,
    options: &Options,
    writer: &mut bcf::Writer,
) -> Result<()> {
    info!(logger, "Processing {}", &path);

    let mut reader = bcf::Reader::from_path(path)
        .chain_err(|| format!("Could not open input file {} for reading", path))?;
    if options.io_threads > 0 {
        reader
            .set_threads(options.io_threads as usize)
            .chain_err(|| "Could not set number of threads")?;
    }

    // If pedigree is given then only founders are considered per-population.
    let considered_samples = if let Some(pedigree) = pedigree {
        let reader_samples = reader
            .header()
            .samples()
            .iter()
            .map(|s| String::from_utf8(s.to_vec()).expect("Could not decode sample name as UTF-8"))
            .collect::<HashSet<_>>();
        pedigree
            .individuals
            .iter()
            .filter(|i| i.father == "0" && i.mother == "0" && reader_samples.contains(&i.name))
            .map(|i| i.name.clone())
            .collect::<HashSet<_>>()
    } else {
        reader
            .header()
            .samples()
            .iter()
            .map(|s| String::from_utf8(s.to_vec()).expect("Could not decode sample name as UTF-8"))
            .collect::<HashSet<_>>()
    };

    let mut i: usize = 0;
    loop {
        let mut record = reader.empty_record();
        match reader.read(&mut record) {
            Ok(_) => record.unpack(),
            Err(bcf::ReadError::NoMoreRecord) => break,
            _ => bail!("Error reading BCF record"),
        }

        // Do not consider SVs.
        if String::from_utf8(
            record
                .info(b"VT")
                .string()
                .unwrap_or_default()
                .unwrap_or_default()
                .pop()
                .unwrap_or_default()
                .to_vec(),
        ).chain_err(|| "Could not decode string")? == "SV"
        {
            continue;
        }

        // Create overall counter.
        let num_alleles = (record.allele_count() - 1) as usize;
        let mut all_stats = GroupStats::new(&"ALL", num_alleles);

        // Create counters for populations.
        let mut pop_stats = pop_map.map(|_| {
            // Initialize per-population counters.
            let mut pop_stats: HashMap<String, GroupStats> = HashMap::new();
            for pop in populations.unwrap() {
                pop_stats.insert(pop.clone(), GroupStats::new(&pop.clone(), num_alleles));
            }
            pop_stats
        });

        // Count genotypes etc.
        {
            let gts = record
                .genotypes()
                .chain_err(|| "Could not get genotypes from record")?;
            for sample in &considered_samples {
                let sample_id = reader
                    .header()
                    .sample_to_id(sample.as_bytes())
                    .chain_err(|| "Sample not found")?;
                let gt = gts.get(*sample_id as usize);

                all_stats.tally(&gt);
                pop_stats.as_mut().map(|map| {
                    let pop_oth = "Oth".to_string();
                    let pop = pop_map
                        .unwrap()
                        .get(sample)
                        .or(Some(&pop_oth))
                        .unwrap();
                    map.get_mut(pop).map(|stats| stats.tally(&gt))
                });
            }
        }

        // Make record suitable for writing out.
        writer.translate(&mut record);
        writer.subset(&mut record);
        trace!(
            logger,
            "{:?} {:?} {:?}",
            &all_stats,
            &all_stats.ac(),
            &all_stats.af()
        );

        // Store count and frequency information in the record.
        record
            .push_info_integer(b"AC", &all_stats.ac())
            .chain_err(|| "Could not write INFO/AC")?;
        record
            .push_info_float(b"AF", &all_stats.af())
            .chain_err(|| "Could not write INFO/AF")?;
        record
            .push_info_integer(b"AN", &[all_stats.an])
            .chain_err(|| "Could not write INFO/AN")?;
        record
            .push_info_integer(b"Hemi", &all_stats.hemi)
            .chain_err(|| "Could not write INFO/Hemi")?;
        record
            .push_info_integer(b"Het", &all_stats.het)
            .chain_err(|| "Could not write INFO/Het")?;
        record
            .push_info_integer(b"Hom", &all_stats.hom)
            .chain_err(|| "Could not write INFO/Hom")?;

        if let Some(pops) = populations {
            let tmp_map = pop_stats.as_ref().unwrap();
            for pop in *pops {
                let stats = tmp_map
                    .get(pop)
                    .expect(&format!("Could not find stats for population {}", &pop));
                record
                    .push_info_integer(format!("{}_AC", &pop).as_bytes(), &stats.ac())
                    .chain_err(|| format!("Could not write INFO/{}_AC", &pop))?;
                record
                    .push_info_float(format!("{}_AF", &pop).as_bytes(), &stats.af())
                    .chain_err(|| format!("Could not write INFO/{}_AF", &pop))?;
                record
                    .push_info_integer(format!("{}_AN", &pop).as_bytes(), &[stats.an])
                    .chain_err(|| format!("Could not write INFO/{}_AN", &pop))?;
                record
                    .push_info_integer(format!("{}_Hemi", &pop).as_bytes(), &stats.hemi)
                    .chain_err(|| format!("Could not write INFO/{}_Hemi", &pop))?;
                record
                    .push_info_integer(format!("{}_Het", &pop).as_bytes(), &stats.het)
                    .chain_err(|| format!("Could not write INFO/{}_Het", &pop))?;
                record
                    .push_info_integer(format!("{}_Hom", &pop).as_bytes(), &stats.hom)
                    .chain_err(|| format!("Could not write INFO/{}_Hom", &pop))?;
            }
            let popmax = pops
                .iter()
                .max_by_key(|p| OrderedFloat(tmp_map.get(*p).unwrap().af().iter().sum::<f32>()))
                .expect("One must be largest");
            let stats = tmp_map
                .get(popmax)
                .expect(&format!("Could not find stats for population {}", &popmax));
            record
                .push_info_string(b"POPMAX", &[popmax.as_bytes()])
                .chain_err(|| "Could not write INFO/POPMAX_AC")?;
            record
                .push_info_integer(b"POPMAX_AC", &stats.ac())
                .chain_err(|| "Could not write INFO/POPMAX_AC")?;
            record
                .push_info_float(b"POPMAX_AF", &stats.af())
                .chain_err(|| "Could not write INFO/POPMAX_AF")?;
            record
                .push_info_integer(b"POPMAX_AN", &[stats.an])
                .chain_err(|| "Could not write INFO/POPMAX_AN")?;
            record
                .push_info_integer(b"POPMAX_Hemi", &stats.hemi)
                .chain_err(|| "Could not write INFO/POPMAX_Hemi")?;
            record
                .push_info_integer(b"POPMAX_Het", &stats.het)
                .chain_err(|| "Could not write INFO/POPMAX_Het")?;
            record
                .push_info_integer(b"POPMAX_Hom", &stats.hom)
                .chain_err(|| "Could not write INFO/POPMAX_Hom")?;
        }

        // Actually write out the record
        writer
            .write(&record)
            .chain_err(|| "Problem writing record to output file")?;
        i += 1;
        if i % 1000 == 0 {
            debug!(logger, "Processed {} records -- at {}", i, record.pos());
        }
    }

    Ok(())
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
    {
        let mut writer = {
            let reader = bcf::Reader::from_path(&options.input[0]).chain_err(|| {
                format!(
                    "Could not open input file {} for reading",
                    &options.input[0]
                )
            })?;
            build_bcf_writer(
                &options.output,
                &populations.as_ref(),
                &contigs.as_ref(),
                &reader,
            )?
        };

        info!(logger, "Starting processing");
        for path in &options.input {
            process_input(
                &logger,
                path,
                &pedigree.as_ref(),
                &pop_map.as_ref(),
                &populations.as_ref(),
                &options,
                &mut writer,
            )?;
        }
    }

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
