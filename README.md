# Summit

Summit is an R-based workflow for annotating SNPs against regulatory genomics resources and transcription factor motif databases.

It takes variants in `chr:pos:ref:alt` format, maps them to rsIDs, queries a chromosome-split tabix-indexed regulatory features bundle, optionally scans allele-specific sequences with FIMO, and writes structured results to a run-specific output directory.

## Main Features

- regulatory overlap queries against a merged tabix-indexed bundle built from:
  - ReMap 2022
  - ChIP-Atlas
  - DNase HINT footprints
- optional FIMO motif scanning with:
  - JASPAR
  - CIS-BP
- allele-specific ref/alt motif comparison tables
- optional overlap against a user-supplied custom tabular dataset
- reproducible config-driven runs

## Inputs

Summit accepts SNPs in:

- `--variants`
- `--variants-file`
- `--config`

Accepted format:

```text
chr1:45439711:G:T
chr1:45480226:C:T
```

Example test variants:

```text
chr1:45439711:G:T
chr1:45480226:C:T
chr1:45506087:T:C
chr1:45540036:G:T
```

## Requirements

To run Summit, users need:

- R with the required CRAN/Bioconductor packages used by `scripts/run_summit.R`
- MEME Suite `fimo` on `PATH` if `run_fimo=true`
- the `merged_regulatory_features` bundle installed under `data/merged_regulatory_features_bgz`
- local motif resource files for FIMO:
  - JASPAR MEME
  - JASPAR TRANSFAC
  - CIS-BP MEME
  - CIS-BP TF information

The bundled regulatory features and built-in references are expected to use hg38 coordinates.

## Download the Regulatory Features Bundle

From the repository root:

```bash
bash scripts/download_database.sh --dest-root .
```

This installs the bundle to:

```text
data/merged_regulatory_features_bgz
```

## Quick Start

Run Summit with the example merged-regulatory-features config:

```bash
Rscript scripts/run_summit.R --config=config/example_merged_regulatory_features.yml
```

Run with explicit variants:

```bash
Rscript scripts/run_summit.R \
  --variants="chr1:45439711:G:T,chr1:45480226:C:T,chr1:45506087:T:C" \
  --merged-regulatory-features-bgz-dir=data/merged_regulatory_features_bgz \
  --output-dir=output/test_run \
  --run-fimo=false
```

## Example Configs

- `config/example_merged_regulatory_features.yml`
- `config/example_no_custom_dataset.yml`
- `config/example_with_custom_dataset.yml`
- `config/example_variants.txt`

## Important Parameters

- `--output-dir`
- `--run-id`
- `--variants`
- `--variants-file`
- `--config`
- `--run-fimo`
- `--force-fimo`
- `--flank-bp`
- `--merged-regulatory-features-bgz-dir`
- `--merged-regulatory-features-prefix`

For custom datasets:

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

## Outputs

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
- optional custom dataset overlap table
- motif metadata tables
- allele-specific sequence tables when `run_fimo=true`
- cached FIMO result files
- FIMO hit tables
- ref/alt motif comparison tables
- prioritized motif/regulatory support tables
- run summary table
- diagnostics table
- session info
- pipeline log

## Notes

- Sequence extraction only happens when `run_fimo=true`.
- If `run_fimo=false`, Summit skips allele-sequence generation and FIMO output.
- The `merged_regulatory_features` bundle preserves source provenance through `source_dataset`.

## Repository Files

- [README.md](/Volumes/One_Touch/GWAS/MUTYH/README.md)
- [scripts/run_summit.R](/Volumes/One_Touch/GWAS/MUTYH/scripts/run_summit.R)
- [scripts/download_database.sh](/Volumes/One_Touch/GWAS/MUTYH/scripts/download_database.sh)
- [config/example_merged_regulatory_features.yml](/Volumes/One_Touch/GWAS/MUTYH/config/example_merged_regulatory_features.yml)
- [config/example_no_custom_dataset.yml](/Volumes/One_Touch/GWAS/MUTYH/config/example_no_custom_dataset.yml)
- [config/example_with_custom_dataset.yml](/Volumes/One_Touch/GWAS/MUTYH/config/example_with_custom_dataset.yml)
- [config/example_variants.txt](/Volumes/One_Touch/GWAS/MUTYH/config/example_variants.txt)
- [license.txt](/Volumes/One_Touch/GWAS/MUTYH/docs/license.txt)
