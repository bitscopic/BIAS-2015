# Changelog

All notable changes to BIAS-2015 will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **Derived missense-O/E × MOI AF cutoffs for BA1/BS1/PM2.** BIAS-2015 now consults a data-derived, constraint-stratified AF cutoff table (`data/af_cutoffs/mis_oe_upper_binned_cutoffs.tsv`) as the second tier of the BA1/BS1/PM2 cascade. Cutoffs are binned by gnomAD v4.1.1 `oe_mis_upper` and stratified by mode of inheritance (AD vs AR).
- **Gene → MOI lookup.** New committed lookup `data/af_cutoffs/gene_to_moi.tsv` resolves per-gene mode of inheritance using HPO → GenCC → ClinGen as ordered fallback sources (expanding the previous ClinGen-only resolution). Loaded by `src/bias_2015/mis_oe_moi_lookup.py`.
- **New `mis_oe_moi_lookup` module.** Provides `load_gene_to_moi`, `load_mis_oe_moi_af_cutoffs`, `resolve_moi`, and `get_mis_oe_moi_cutoff` for the derived-cutoff branch of BA1/BS1/PM2.
- **Calibration pipeline for the derived cutoffs.** New scripts under `src/scripts/calibration/` (`05`–`09`) that filter ClinVar/Nirvana, join MOI and gnomAD constraints, derive constraint-stratified cutoffs via KDE + Tavtigian LR thresholds with cluster bootstrap CIs, and produce cross-validation plots.

### Changed

- **BA1/BS1/PM2 evaluation cascade.** Population-AF codes now resolve in tiered order: VCEP rule (if present) → derived missense-O/E × MOI cutoff → data-derived flat fallback in `constants.py`. The flat fallback replaces the historical ACMG 2015 defaults (BA1 5%, BS1 0.1%, PM2 0.1%), which were empirically too permissive during ClinVar validation.
- **BS1 fires at three strengths (Strong/Moderate/Supporting).** Matches the stratification of the derived cutoff table.

## [3.0.0] - 2026-03-05

### Added

- **VEP annotator support.** BIAS-2015 now supports Ensembl VEP 115.0 as an alternative to Nirvana 3.18.1 for variant annotation. A new `--annotator` CLI flag selects between `nirvana` (default) and `vep`. New module `src/bias_2015/extract_from_vep_json.py` handles VEP JSON parsing and variant classification. Full setup instructions are provided in `doc/vep_setup.md`.

- **AlphaMissense integration.** Missense variant pathogenicity scores from AlphaMissense are now used in PP3 and BP4 classifiers. PP3 thresholds: supporting 0.792, moderate 0.906, strong 0.972, very strong 0.99. BP4 thresholds: pathogenic cutoff 0.792, tiered benign cutoffs 0.17, 0.10, 0.07. Thresholds follow Bergquist et al. *Genet Med* 2025. Available for both Nirvana (via custom supplementary annotation) and VEP (via plugin).

- **ClinVar submitter counts for PS4.** New `generate_clinvar_submitter_counts.py` preprocessing script produces per-variant submitter counts from ClinVar XML, used by the PS4 classifier as a fallback when GWAS data is unavailable. Thresholds: supporting (2), moderate (4), strong (6).

- **Dedicated variant data model.** New `src/bias_2015/variant_data.py` provides the `VariantData` class — a shared data structure used by both Nirvana and VEP extraction modules with `to_tsv()` and `to_json()` output methods.

- **Centralized gene data loader.** New `src/bias_2015/gene_data_loader.py` loads ClinGen Gene-Disease Validity and gnomAD gene constraint data independently of the annotator, providing a single source of truth for gene-level information.

- **Shared output utilities.** New `src/scripts/bias_output_utils.py` provides common functions for reading and comparing BIAS output files.

- **Annotator-specific test data.** Test directory now includes separate Nirvana and VEP input/output files (`bias-2015_test_file.nirvana.json`, `bias-2015_test_file.vep.json`, and corresponding `.bias_output.tsv` files) for 100 randomly selected eRepo variants.

- **Updated eRepo VCF.** New `erepo_03.02.2026.vcf` generated from the latest eRepo dataset with improved indel coordinate conversion.

- **VEP setup documentation.** New `doc/vep_setup.md` provides complete instructions for installing VEP via Docker, downloading required annotation sources (gnomAD, ClinVar, REVEL, AlphaMissense), and running VEP preprocessing.


### Changed

- **Dual-annotator architecture.** The `classify_variants()` function has been moved from the main script into annotator-specific modules (`extract_from_nirvana_json.py` and `extract_from_vep_json.py`). The main `bias_2015.py` entry point dispatches to the appropriate module based on the `--annotator` flag.

- **PS4 TopMed fallback replaced with ClinVar submitter counts.** The PS4 classifier no longer uses TopMed allele/homozygous counts as a fallback. Instead, it uses ClinVar pathogenic submitter counts when GWAS data is unavailable.

- **ClinVar entry selection improvements.** New `normalize_alleles()` and `alleles_match()` functions handle VCF padding vs Nirvana-style allele normalization. ClinVar entries without significance assertions are now filtered out.

- **Preprocessing unified across annotators.** `preprocessing.py` now supports both annotators via `--annotator nirvana|vep|all` and includes a `skip_if_exists()` helper to avoid redundant downloads and processing steps. VEP annotation of ClinVar is performed via Docker during preprocessing.

- **Updated preprocessing data (March 1, 2026).** All BIAS data files have been regenerated using ClinVar February 2026 and updated supplementary annotations. Data files are hosted on AWS S3 at `s3://bias-2015/v3.0.0_datasets/2026.03.01/`.

- **Updated Nirvana supplementary annotations.** Updated ClinVar (February 2026) and new AlphaMissense NSA files are hosted on S3 at `s3://bias-2015/nirvana_data/` to replace stale files bundled with Nirvana 3.18.1.

- **Test data reorganized.** Old monolithic test files have been replaced with annotator-specific test inputs and expected outputs, enabling deterministic validation for both Nirvana and VEP pipelines.

- **README updated.** The README now documents multi-annotator setup, S3 data downloads, VEP as an alternative annotator, and updated preprocessing instructions.

### Fixed

- **ClinVar indel parsing.** Fixed incorrect coordinate handling for ClinVar insertions and deletions via new allele normalization, improving variant matching accuracy in PS1/PM5 and other ClinVar-dependent classifiers.

- **Transcript assignment bug.** Resolved an issue where the wrong transcript could be selected for variants with multiple transcript annotations. Transcript ranking now uses a weighted system prioritizing canonical status (+3), protein-coding bioType (+2), and coding sequence information (+2) over database source (+1) and HGNC mapping (+1).

- **ABSplice off by one bug.** Fixed ABSplice preprocessing which was causing an off by one coordinate mapping issue. 

### Removed

- **`variant_interpretation_dataset_loader.py` deleted.** The monolithic dataset loader has been replaced by the modular `bias_dataset_loader.py`, `gene_data_loader.py`, and `variant_data.py` modules.

- **Old monolithic test files.** Legacy test data files (`bias-2015-expected_test_output.tsv`, `bias-2015_test_file.json`) that did not distinguish between annotators have been removed in favor of annotator-specific test data.

### Upgrade Notes

Users upgrading from v2.1.1 should note:

1. **New data files required.** v2.x.x preprocessing data files are not compatible with v3.0.0. Download new data from `s3://bias-2015/v3.0.0_datasets/2026.03.01/` or regenerate via `preprocessing.py`.
2. **Updated `required_paths.json` format.** The configuration file includes additional paths for AlphaMissense and ClinVar submitter count data.
3. **Updated Nirvana supplementary annotations.** Install updated ClinVar and AlphaMissense NSA files from `s3://bias-2015/nirvana_data/`.
4. **New CLI flag.** The `--annotator` flag defaults to `nirvana` — existing Nirvana workflows are unaffected without code changes.

## [2.1.1] - 2025-11-01

Initial public release accompanying Eisenhart et al. "Automating ACMG variant classifications
with BIAS-2015 v2.1.1: algorithm analysis and benchmark against the FDA-approved eRepo dataset."
*Genome Medicine* 17, 148 (2025).
