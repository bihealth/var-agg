//! I/O for pedigree files.

use std::fs::File;

use csv;

use super::errors::*;

/// Sex of an individual
#[derive(Clone, Debug, Copy, PartialEq)]
pub enum Sex {
    Male,
    Female,
    Unknown,
}

impl Sex {
    /// Return `Sex` from PED file integer.
    pub fn from_int(i: i32) -> Self {
        match i {
            1 => Sex::Male,
            2 => Sex::Female,
            _ => Sex::Unknown,
        }
    }

    /// Return integer representation for PED file.
    #[allow(dead_code)]
    pub fn to_int(&self) -> i32 {
        match self {
            Sex::Male => 1,
            Sex::Female => 2,
            Sex::Unknown => 0,
        }
    }
}

/// Affection status of an individual.
#[derive(Clone, Debug, Copy, PartialEq)]
pub enum DiseaseStatus {
    Affected,
    Unaffected,
    Unknown,
}

impl DiseaseStatus {
    /// Return `DiseaseStatus` from PED file integer.
    pub fn from_int(i: i32) -> Self {
        match i {
            1 => DiseaseStatus::Unaffected,
            2 => DiseaseStatus::Affected,
            _ => DiseaseStatus::Unknown,
        }
    }

    /// Return integer representation for PED file.
    #[allow(dead_code)]
    pub fn to_int(&self) -> i32 {
        match self {
            DiseaseStatus::Unaffected => 1,
            DiseaseStatus::Affected => 2,
            DiseaseStatus::Unknown => 0,
        }
    }
}

/// Store information for one person.
#[derive(Clone, Debug)]
pub struct Individual {
    /// ID of the family
    pub family: String,
    /// ID of the individual
    pub name: String,
    /// ID of the father
    pub father: String,
    /// ID of the mother
    pub mother: String,
    /// The sex of the individual
    pub sex: Sex,
    /// The disease status of the individual
    pub disease_status: DiseaseStatus,
    /// Any further fields from the PED files, tokenized at TAB characters
    pub extra_fields: Vec<String>,
}

/// Store information for one pedigree.
#[derive(Clone, Debug)]
pub struct Pedigree {
    /// The individuals in the pedigree
    pub individuals: Vec<Individual>,
    /// The headings of the core fields (up to disease status)
    pub core_headings: Vec<String>,
    /// The headings of the extra files (starting after disease status)
    pub extra_headings: Vec<String>,
}

impl Pedigree {
    /// Load `Pedigree` from path to PED file.
    pub fn from_path(path: &str) -> Result<Pedigree> {
        // Create CSV reader and check columns.
        let file =
            File::open(&path).chain_err(|| format!("Could not open input file: {}", &path))?;
        let mut reader = csv::ReaderBuilder::new()
            .flexible(true)
            .delimiter(b'\t')
            .has_headers(true)
            .from_reader(file);
        if reader
            .headers()
            .chain_err(|| "Could not access PED headers")?
            .len()
            < 6
        {
            bail!("PED file has less than 6 columns!");
        }

        // Extract headings.
        let core_headings = reader
            .headers()
            .chain_err(|| "Could not access headers")?
            .iter()
            .take(6)
            .map(|s| s.to_string())
            .collect::<Vec<_>>();
        let extra_headings = reader
            .headers()
            .chain_err(|| "Could not access headers")?
            .iter()
            .skip(6)
            .map(|s| s.to_string())
            .collect::<Vec<_>>();

        // Load data
        let mut individuals = Vec::new();
        for record in reader.records() {
            let record = record.chain_err(|| "Problem reading PED record from file")?;
            if record.len() < 6 {
                bail!("PED record had less than 6 columns");
            }
            individuals.push(Individual {
                family: record.get(0).unwrap().to_string(),
                name: record.get(1).unwrap().to_string(),
                father: record.get(2).unwrap().to_string(),
                mother: record.get(3).unwrap().to_string(),
                sex: Sex::from_int(
                    record
                        .get(4)
                        .unwrap()
                        .parse::<i32>()
                        .expect("Could not parse column \"sex\""),
                ),
                disease_status: DiseaseStatus::from_int(
                    record
                        .get(5)
                        .unwrap()
                        .parse::<i32>()
                        .expect("Could not parse column \"disease_status\""),
                ),
                extra_fields: record
                    .iter()
                    .skip(6)
                    .map(|s| s.to_string())
                    .collect::<Vec<_>>(),
            });
        }

        Ok(Pedigree {
            individuals,
            core_headings,
            extra_headings,
        })
    }
}
