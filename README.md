# Summit

## Overview

This repository contains Summit, an R-based workflow for annotating SNPs against regulatory genomics resources and motif databases.

The main entry point is:

- `scripts/run_pipeline.R`

This script is the reproducible Summit workflow version of the earlier exploratory workflow in:

- `scripts/MUTYH_workflow_clean_031226.R`

Summit accepts user-supplied SNPs in `chr:pos:ref:alt` format, maps them to rsIDs, queries a merged tabix-indexed regulatory features bundle, optionally searches a custom dataset, optionally runs FIMO motif scans, and writes structured outputs to a per-run directory.

## Current Built-In Regulatory Query Layer

Summit expects a normalized merged regulatory features bundle split by chromosome and stored as:

- `merged_regulatory_features_chr1.tsv.bgz`
- `merged_regulatory_features_chr1.tsv.bgz.tbi`
- `merged_regulatory_features_chr2.tsv.bgz`
- `merged_regulatory_features_chr2.tsv.bgz.tbi`
- ...

These files are queried through `Rsamtools`/`tabix` and replace the older direct runtime reads of:

- ReMap 2022
- ChIP-Atlas
- DNase HINT 16 bp footprints
- DNase HINT 20 bp footprints

Dataset identity is preserved in the merged bundle with `source_dataset`, so Summit still emits separate overlap tables for:

- ReMap
- ChIP-Atlas
- DNase HINT16
- DNase HINT20

Additional built-in resources still include:

- dbSNP155 for rsID lookup
- hg38 reference genome for sequence extraction
- JASPAR motif database for FIMO
- CIS-BP motif database for FIMO
- a FIMO background model

## Main Inputs

Summit accepts SNPs in these ways:

- inline with `--variants`
- from a text file with `--variants-file`
- from a simple config file with `--config`

Accepted SNP format:

```text
chr1:45439711:G:T
chr1:45480226:C:T
```

Example SNPs you can use for quick testing:

```text
chr1:45439711:G:T
chr1:45480226:C:T
chr1:45506087:T:C
chr1:45540036:G:T
```

## Config Examples

Example config files in `config/`:

- `config/example_no_custom_dataset.yml`
- `config/example_with_custom_dataset.yml`
- `config/example_merged_regulatory_features.yml`
- `config/example_variants.txt`

The merged-bundle example is the best starting point for Summit:

```bash
Rscript scripts/run_pipeline.R --config=config/example_merged_regulatory_features.yml
```

By default, Summit looks for the bundled interval and motif resource files under the repo’s `data/` directory. If `run_fimo=true`, it will also try to auto-detect the MEME Suite `fimo` executable from your `PATH`; set `meme_bin` in config only if auto-detection is not sufficient.

Users can download the bundled `merged_regulatory_features` resource from Zenodo with:

```bash
bash scripts/download_database.sh --dest-root .
```

## Important Parameters

Common parameters:

- `--output-dir`
- `--run-id`
- `--variants`
- `--variants-file`
- `--config`
- `--run-fimo`
- `--force-fimo`
- `--flank-bp`

Merged regulatory features bundle parameters:

- `--merged-regulatory-features-bgz-dir`
- `--merged-regulatory-features-prefix`

Custom dataset parameters:

- `--custom-dataset-name`
- `--custom-dataset-file`
- `--custom-dataset-files`
- `--custom-dataset-glob`
- `--custom-dataset-chr-col`
- `--custom-dataset-start-col`
- `--custom-dataset-end-col`
- `--custom-dataset-delim`
- `--custom-dataset-has-header`
- `--custom-dataset-zero-based`

## Example Usage

Run with the merged regulatory features bundle from config:

```bash
Rscript scripts/run_pipeline.R --config=config/example_merged_regulatory_features.yml
```

Run with explicit SNPs and a merged regulatory features bundle directory:

```bash
Rscript scripts/run_pipeline.R \
  --variants="chr1:45439711:G:T,chr1:45480226:C:T,chr1:45506087:T:C" \
  --merged-regulatory-features-bgz-dir=/path/to/merged_regulatory_features_bgz \
  --output-dir=output/test_run \
  --run-fimo=false
```

Run with a custom dataset as well:

```bash
Rscript scripts/run_pipeline.R \
  --config=config/example_merged_regulatory_features.yml \
  --custom-dataset-name=my_peaks \
  --custom-dataset-file=/path/to/custom_dataset.tsv.gz \
  --custom-dataset-chr-col=chrom \
  --custom-dataset-start-col=start \
  --custom-dataset-end-col=end
```

## Custom Dataset Support

Users can add one extra tabular dataset, or a collection of tabular files, to search against the same SNPs.

Supported modes:

- one file: `--custom-dataset-file=/path/to/file.tsv.gz`
- many files: `--custom-dataset-files=file1,file2,file3`
- glob pattern: `--custom-dataset-glob='/path/to/chr*.tsv.gz'`

Required metadata:

- chromosome column, unless chromosome can be inferred from filenames
- start column
- end column

Optional metadata:

- delimiter: `auto`, `tsv`, or `csv`
- header present or absent
- whether starts are zero-based

## Output Structure

Each run creates:

```text
<output_dir>/
  config/
  diagnostics/
  logs/
  results/
    tables/
    rds/
```

Typical outputs include:

- input SNP table
- per-dataset overlap tables
- overlap table for the custom dataset, if used
- motif metadata tables
- allele-specific sequence tables when `run_fimo=true`
- cached FIMO result files
- FIMO hit tables
- run summary table
- diagnostics table
- session info
- pipeline log

## Sequence and FIMO Behavior

Sequence extraction now only happens when `run_fimo=true`.

So:

- `run_fimo=false`: no flank sequence generation and no FIMO output
- `run_fimo=true`: build allele-specific flank sequences and run FIMO

## Coordinate Assumption

The internal regulatory interval bundle and built-in reference resources are expected to use hg38 coordinates. Input SNPs should therefore also be provided on hg38 coordinates.

## Interval Bundle Build Scripts

The repository also contains preprocessing scripts for building the normalized merged regulatory features bundle, including:

- `scripts/prepare_chip_atlas_metadata_table_ultrafast.R`
- `scripts/prepare_chip_atlas_normalized_by_chr_ultrafast.R`
- `scripts/prepare_remap_normalized_by_chr_ultrafast.R`
- `scripts/prepare_dnase_hint_normalized_by_chr_ultrafast.R`
- `scripts/merge_normalized_intervals_by_chr.R`
- `scripts/index_merged_intervals_bgzip_tabix.R`
- `scripts/spot_check_merged_normalized_by_chr.R`

There is also a wrapper for the speed-first HPC workflow:

- `scripts/run_all_normalization_ultrafast.R`

## Regulatory Features Bundle Design

The normalized interval bundle schema/design notes are documented in:

- [docs/normalized_interval_dataset_spec.md](/Volumes/One_Touch/GWAS/MUTYH/docs/normalized_interval_dataset_spec.md)

That spec captures the shared schema used across ReMap, ChIP-Atlas, and DNase HINT, including dedicated fields such as:

- `experiment_type`
- `cell_type`
- `tissue_type`
- `gene`
- `motif`

## Repository Files

- [scripts/run_pipeline.R](/Volumes/One_Touch/GWAS/MUTYH/scripts/run_pipeline.R)
- [scripts/download_database.sh](/Volumes/One_Touch/GWAS/MUTYH/scripts/download_database.sh)
- [scripts/run_all_normalization_ultrafast.R](/Volumes/One_Touch/GWAS/MUTYH/scripts/run_all_normalization_ultrafast.R)
- [scripts/merge_normalized_intervals_by_chr.R](/Volumes/One_Touch/GWAS/MUTYH/scripts/merge_normalized_intervals_by_chr.R)
- [scripts/index_merged_intervals_bgzip_tabix.R](/Volumes/One_Touch/GWAS/MUTYH/scripts/index_merged_intervals_bgzip_tabix.R)
- [config/example_merged_regulatory_features.yml](/Volumes/One_Touch/GWAS/MUTYH/config/example_merged_regulatory_features.yml)
- [config/example_no_custom_dataset.yml](/Volumes/One_Touch/GWAS/MUTYH/config/example_no_custom_dataset.yml)
- [config/example_with_custom_dataset.yml](/Volumes/One_Touch/GWAS/MUTYH/config/example_with_custom_dataset.yml)
- [config/example_variants.txt](/Volumes/One_Touch/GWAS/MUTYH/config/example_variants.txt)
- [docs/normalized_interval_dataset_spec.md](/Volumes/One_Touch/GWAS/MUTYH/docs/normalized_interval_dataset_spec.md)
- [license.txt](/Volumes/One_Touch/GWAS/MUTYH/docs/license.txt)
- [notes.md](/Volumes/One_Touch/GWAS/MUTYH/notes.md)
