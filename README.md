# Var-Agg

A simple helper for aggregating multi-sample VCF files into "site VCF" files,
e.g., for aggregating the thousand genomes per-sample genotype files or in-house
cohorts.

Note that the frequencies are computed for founders only (both mother and father
are "0" in the pedigree file).  Counts including non-founders are in the implicit
"ALL" cohort.
