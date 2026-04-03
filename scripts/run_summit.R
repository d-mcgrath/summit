###############################################################################
# Reproducible regulatory annotation pipeline
###############################################################################

required_packages <- c(
  "BSgenome",
  "BSgenome.Hsapiens.UCSC.hg38",
  "SNPlocs.Hsapiens.dbSNP155.GRCh38",
  "GenomicRanges",
  "Rsamtools",
  "Biostrings",
  "memes",
  "universalmotif",
  "tidyverse"
)

bioc_packages <- c(
  "BSgenome",
  "BSgenome.Hsapiens.UCSC.hg38",
  "SNPlocs.Hsapiens.dbSNP155.GRCh38",
  "GenomicRanges",
  "Rsamtools",
  "Biostrings"
)

cran_packages <- c(
  "memes",
  "universalmotif",
  "tidyverse"
)

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing_packages) > 0) {
  missing_bioc <- intersect(missing_packages, bioc_packages)
  missing_cran <- intersect(missing_packages, cran_packages)

  lines <- c(
    "Summit cannot start because required R packages are missing.",
    "",
    paste0("Missing packages: ", paste(missing_packages, collapse = ", "))
  )

  if (length(missing_bioc) > 0) {
    lines <- c(
      lines,
      "",
      "Bioconductor packages to install:",
      paste0("  ", paste(missing_bioc, collapse = ", ")),
      paste0(
        "  Example: BiocManager::install(c(",
        paste(sprintf('"%s"', missing_bioc), collapse = ", "),
        "))"
      )
    )
  }

  if (length(missing_cran) > 0) {
    lines <- c(
      lines,
      "",
      "CRAN packages to install:",
      paste0("  ", paste(missing_cran, collapse = ", ")),
      paste0(
        "  Example: install.packages(c(",
        paste(sprintf('"%s"', missing_cran), collapse = ", "),
        "))"
      )
    )
  }

  lines <- c(lines, "", "Install the missing package(s) before running scripts/run_summit.R.")
  writeLines(lines, con = stderr())
  quit(save = "no", status = 1)
}

suppressPackageStartupMessages({
  library(BSgenome)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
  library(GenomicRanges)
  library(Rsamtools)
  library(Biostrings)
  library(memes)
  library(universalmotif)
  library(tidyverse)
})

script_file_arg <- commandArgs(FALSE)
script_file <- script_file_arg[grepl("^--file=", script_file_arg)]
script_path <- if (length(script_file) > 0) {
  normalizePath(sub("^--file=", "", script_file[[1]]), mustWork = FALSE)
} else {
  normalizePath(file.path(getwd(), "scripts", "run_summit.R"), mustWork = FALSE)
}
script_dir <- dirname(script_path)
project_root <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
default_meme_bin <- {
  fimo_path <- Sys.which("fimo")
  if (nzchar(fimo_path)) dirname(fimo_path) else NULL
}

###############################################################################
# Defaults
###############################################################################

default_paths <- list(
  meme_bin = default_meme_bin,
  jaspar_motifs = file.path(project_root, "data", "fimo_motif_data", "JASPAR2024_CORE_non-redundant_pfms_meme.txt"),
  jaspar_transfac = file.path(project_root, "data", "fimo_motif_data", "JASPAR2024_CORE_vertebrates_non-redundant_pfms_transfac.txt"),
  cisbp_meme = file.path(project_root, "data", "fimo_motif_data", "CIS-BP", "cisbp_all_motifs.meme"),
  cisbp_info = file.path(project_root, "data", "fimo_motif_data", "CIS-BP", "TF_Information_all_motifs_plus.txt"),
  fimo_bg = file.path(project_root, "data", "fimo_background_models", "fimo_bg_model_2nd_order.txt"),
  merged_regulatory_features_bgz_dir = file.path(project_root, "data", "merged_regulatory_features_bgz"),
  merged_regulatory_features_prefix = "merged_regulatory_features"
)

if (!file.exists(default_paths$fimo_bg)) {
  default_paths$fimo_bg <- NULL
}

###############################################################################
# CLI parsing
###############################################################################

parse_cli_args <- function(args) {
  parsed <- list()

  for (arg in args) {
    if (!startsWith(arg, "--")) {
      stop("All arguments must use --name=value or --flag syntax.")
    }

    body <- sub("^--", "", arg)

    if (grepl("=", body, fixed = TRUE)) {
      key <- sub("=.*$", "", body)
      value <- sub("^[^=]+=", "", body)
      parsed[[key]] <- value
    } else {
      parsed[[body]] <- TRUE
    }
  }

  parsed
}

cli_args <- parse_cli_args(commandArgs(trailingOnly = TRUE))

parse_scalar_config_value <- function(value) {
  trimmed <- str_trim(value)

  if (trimmed %in% c("", "null", "NULL", "~")) {
    return(NULL)
  }

  if (tolower(trimmed) %in% c("true", "false")) {
    return(tolower(trimmed) == "true")
  }

  if (str_detect(trimmed, "^-?\\d+$")) {
    return(as.integer(trimmed))
  }

  if (
    (str_starts(trimmed, "\"") && str_ends(trimmed, "\"")) ||
    (str_starts(trimmed, "'") && str_ends(trimmed, "'"))
  ) {
    return(str_sub(trimmed, 2, -2))
  }

  trimmed
}

read_simple_config <- function(path) {
  lines <- readLines(path, warn = FALSE)
  config <- list()

  for (line in lines) {
    trimmed <- str_trim(line)

    if (trimmed == "" || str_starts(trimmed, "#")) {
      next
    }

    if (!str_detect(trimmed, ":")) {
      next
    }

    key <- str_trim(str_replace(trimmed, ":.*$", ""))
    value <- str_trim(str_replace(trimmed, "^[^:]+:", ""))
    config[[key]] <- parse_scalar_config_value(value)
  }

  config
}

resolve_input_path <- function(path, base_dirs = character(), must_exist = FALSE, path_type = c("file", "dir")) {
  path_type <- match.arg(path_type)

  if (is.null(path) || identical(path, "")) {
    return(path)
  }

  path <- as.character(path)

  path_exists <- function(candidate) {
    if (path_type == "dir") {
      dir.exists(candidate)
    } else {
      file.exists(candidate)
    }
  }

  is_absolute <- function(x) {
    grepl("^(/|[A-Za-z]:[/\\\\])", x)
  }

  candidates <- if (is_absolute(path)) {
    path
  } else {
    unique(c(path, file.path(base_dirs, path)))
  }

  for (candidate in candidates) {
    if (path_exists(candidate)) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
  }

  if (must_exist) {
    stop("Required ", path_type, " not found: ", path)
  }

  normalizePath(candidates[[length(candidates)]], mustWork = FALSE)
}

config_path_cli <- cli_args[["config"]]
config_path_cli <- if (!is.null(config_path_cli)) {
  resolve_input_path(
    config_path_cli,
    base_dirs = unique(c(getwd(), project_root)),
    must_exist = TRUE,
    path_type = "file"
  )
} else {
  NULL
}
config_args <- if (!is.null(config_path_cli)) read_simple_config(config_path_cli) else list()

get_arg <- function(name, default = NULL, config_aliases = character()) {
  if (!is.null(cli_args[[name]])) {
    return(cli_args[[name]])
  }

  config_keys <- unique(c(gsub("-", "_", name), config_aliases))

  for (key in config_keys) {
    if (!is.null(config_args[[key]])) {
      return(config_args[[key]])
    }
  }

  default
}

coerce_flag <- function(x, default = FALSE) {
  if (is.null(x)) {
    return(default)
  }

  if (isTRUE(x)) {
    return(TRUE)
  }

  value <- tolower(as.character(x))
  value %in% c("1", "true", "t", "yes", "y")
}

coerce_delim <- function(x, default = "auto") {
  if (is.null(x)) {
    return(default)
  }

  value <- tolower(as.character(x))
  if (!value %in% c("auto", "tsv", "csv")) {
    stop("Delimiter must be one of: auto, tsv, csv")
  }

  value
}

###############################################################################
# Helpers
###############################################################################

timestamp_string <- function() {
  format(Sys.time(), "%Y%m%d_%H%M%S")
}

normalize_chr <- function(chr) {
  chr |>
    as.character() |>
    str_replace("^chr", "") |>
    str_to_upper()
}

format_chr_for_ucsc <- function(chr) {
  paste0("chr", normalize_chr(chr))
}

format_chr_numeric <- function(chr) {
  normalized <- normalize_chr(chr)
  suppressWarnings(as.integer(normalized))
}

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

pipeline_log_con <- NULL

log_message <- function(..., .time = TRUE) {
  prefix <- if (.time) paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ") else ""
  line <- paste0(prefix, paste0(..., collapse = ""))
  cat(line, "\n")

  if (!is.null(pipeline_log_con) && isOpen(pipeline_log_con)) {
    writeLines(line, pipeline_log_con)
    flush(pipeline_log_con)
  }
}

safe_write_tsv <- function(tbl, path) {
  readr::write_tsv(tbl, path, na = "")
  invisible(path)
}

safe_save_rds <- function(object, path) {
  saveRDS(object, path)
  invisible(path)
}

parse_variant_tbl <- function(variant_ids) {
  tibble(variant_id = variant_ids) |>
    mutate(
      variant_id = str_trim(variant_id),
      chr = str_replace(variant_id, "^chr", ""),
      chr = str_replace(chr, ":(.*)$", ""),
      pos = str_replace(variant_id, "^chr", ""),
      pos = str_replace(pos, "^[^:]+:(\\d+):.*$", "\\1"),
      ref = str_replace(variant_id, "^chr", ""),
      ref = str_replace(ref, "^[^:]+:\\d+:([^:]+):.*$", "\\1"),
      alt = str_replace(variant_id, "^chr", ""),
      alt = str_replace(alt, "^[^:]+:\\d+:[^:]+:([^:]+)$", "\\1")
    ) |>
    mutate(
      chr = normalize_chr(chr),
      pos = as.integer(pos),
      ref = str_to_upper(ref),
      alt = str_to_upper(alt)
    )
}

validate_variant_tbl <- function(tbl) {
  invalid <- tbl |>
    filter(
      is.na(pos) |
        !str_detect(chr, "^(\\d+|X|Y|M|MT)$") |
        !str_detect(ref, "^[ACGT]+$") |
        !str_detect(alt, "^[ACGT]+$")
    )

  if (nrow(invalid) > 0) {
    stop(
      "Invalid variant IDs detected. Expected format like chr1:12345:A:G. ",
      "Offending values: ",
      paste(invalid$variant_id, collapse = ", ")
    )
  }

  tbl
}

load_variants <- function(variants_arg = NULL, variants_file = NULL) {
  if (!is.null(variants_file)) {
    raw <- readr::read_lines(variants_file)
    variants <- raw[raw != ""]
  } else if (!is.null(variants_arg)) {
    variants <- str_split(variants_arg, ",", simplify = TRUE) |>
      as.character() |>
      str_trim()
    variants <- variants[variants != ""]
  } else {
    stop(
      "No variants were provided. Supply SNPs with --variants, --variants-file, or a config file key such as variants or variants_file. ",
      "Expected format: chr1:45439711:G:T. Built-in interval resources use hg38 coordinates."
    )
  }

  parse_variant_tbl(variants) |>
    validate_variant_tbl() |>
    distinct(variant_id, .keep_all = TRUE)
}

make_snp_granges <- function(tbl, chr_col = "chr", pos_col = "pos") {
  GRanges(
    seqnames = format_chr_for_ucsc(tbl[[chr_col]]),
    ranges = IRanges(
      start = as.numeric(tbl[[pos_col]]),
      end = as.numeric(tbl[[pos_col]])
    )
  )
}

get_rsids <- function(variant_tbl, snp_db = SNPlocs.Hsapiens.dbSNP155.GRCh38::SNPlocs.Hsapiens.dbSNP155.GRCh38) {
  gr <- make_snp_granges(variant_tbl)
  known_snps <- suppressWarnings(snpsByOverlaps(snp_db, gr))
  hits <- suppressWarnings(findOverlaps(gr, known_snps))
  mapping <- selectHits(hits, select = "first")

  variant_tbl |>
    mutate(
      rsid = mcols(known_snps)$RefSNP_id[mapping],
      rsid = if_else(is.na(rsid), NA_character_, paste0("rs", rsid))
    )
}

flank_seq_for_row <- function(chr, pos, flank = 20) {
  genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  chr_ucsc <- format_chr_for_ucsc(chr)
  chr_len <- seqlengths(genome)[chr_ucsc]
  start0 <- max(1, pos - flank)
  end0 <- min(chr_len, pos + flank)

  tibble(
    flank_start = start0,
    flank_end = end0,
    sequence = as.character(getSeq(genome, chr_ucsc, start = start0, end = end0))
  )
}

add_flank_sequences <- function(variant_tbl, flank = 20) {
  flank_tbl <- pmap_dfr(
    list(variant_tbl$chr, variant_tbl$pos),
    ~ flank_seq_for_row(..1, ..2, flank = flank)
  )

  bind_cols(variant_tbl, flank_tbl)
}

build_allele_sequences <- function(variant_tbl, flank = 20) {
  expected_width <- (flank * 2) + 1

  variant_tbl |>
    add_flank_sequences(flank = flank) |>
    mutate(sequence_width = nchar(sequence)) |>
    filter(sequence_width == expected_width) |>
    select(-sequence_width) |>
    uncount(2) |>
    group_by(variant_id) |>
    mutate(
      allele_index = row_number(),
      allele = if_else(allele_index == 1, "ref", "alt"),
      allele_base = if_else(allele == "ref", ref, alt),
      updated_sequence = {
        tmp <- sequence
        str_sub(tmp, flank + 1, flank + 1) <- allele_base
        tmp
      },
      center_base_original = str_sub(sequence, flank + 1, flank + 1),
      center_base_updated = str_sub(updated_sequence, flank + 1, flank + 1)
    ) |>
    ungroup()
}

build_seq_tbl <- function(rsSeqs, seq_metadata = NULL) {
  seq_tbl <- rsSeqs |>
    imap_dfr(~ tibble(seqnames = .y, sequence = as.character(.x)))

  if (is.null(seq_metadata)) {
    return(seq_tbl)
  }

  seq_tbl |>
    left_join(seq_metadata, by = "seqnames")
}

build_motif_info <- function(cisbp_info_filepath, jaspar_transfac_filepath) {
  cisbp_info_raw <- read_tsv(
    cisbp_info_filepath,
    show_col_types = FALSE,
    name_repair = "minimal"
  )

  cisbp_info <- cisbp_info_raw |>
    select(TF_ID, Motif_ID, TF_Name, Family_Name, DBDs) |>
    rename_with(str_to_lower) |>
    group_by(tf_id, tf_name, family_name, dbds) |>
    summarize(motif_ids = paste0(motif_id, collapse = ";"), .groups = "drop")

  jaspar_info <- tibble(line = readLines(jaspar_transfac_filepath)) |>
    filter(str_detect(line, "^AC ") | str_detect(line, "^CC ") | str_detect(line, "^ID ")) |>
    separate(line, into = c("type", "value"), sep = " ", extra = "merge")

  jaspar_info <- jaspar_info |>
    mutate(
      value = str_replace(value, "::", ";"),
      type = str_to_lower(type),
      type = case_when(
        str_detect(value, ":") ~ str_replace(value, "(.*):.*", "\\1"),
        type == "ac" ~ "tf_id",
        type == "id" ~ "tf_name",
        TRUE ~ type
      )
    ) |>
    pivot_wider(
      names_from = type,
      values_from = value,
      values_fn = \(value) paste(unique(value), collapse = ";")
    ) |>
    mutate(
      tf_name = str_replace(tf_name, ";", "::"),
      tf_family = str_replace(tf_family, "tf_family:", ""),
      tf_class = str_replace(tf_class, "tf_class:", ""),
      tax_group = str_replace(tax_group, "tax_group:", ""),
      pubmed_ids = str_replace(pubmed_ids, "pubmed_ids:", ""),
      uniprot_ids = str_replace(uniprot_ids, "uniprot_ids:", ""),
      data_type = str_replace(data_type, "data_type:", "")
    )

  motif_info <- jaspar_info |>
    mutate(db = "jaspar") |>
    bind_rows(
      cisbp_info |>
        rename(tf_family = family_name) |>
        select(-dbds) |>
        mutate(db = "cis-bp")
    )

  motif_info_long <- jaspar_info |>
    mutate(db = "jaspar", motif_id = tf_id) |>
    bind_rows(
      cisbp_info |>
        rename(tf_family = family_name) |>
        select(-dbds) |>
        mutate(db = "cis-bp") |>
        rename(motif_id = motif_ids) |>
        separate_rows(motif_id, sep = ";")
    )

  list(
    cisbp_info = cisbp_info,
    jaspar_info = jaspar_info,
    motif_info = motif_info,
    motif_info_long = motif_info_long
  )
}

combine_fimo_results <- function(fimo_results_cisbp, fimo_results_jaspar, rsSeqs, motif_info_long, seq_metadata = NULL) {
  bind_rows(
    as_tibble(fimo_results_cisbp) |> mutate(database = "CIS-BP"),
    as_tibble(fimo_results_jaspar) |> mutate(database = "JASPAR")
  ) |>
    left_join(build_seq_tbl(rsSeqs, seq_metadata = seq_metadata), by = "seqnames") |>
    left_join(
      motif_info_long |>
        select(tf_id, tf_name, tf_family, db, motif_id) |>
        group_by(motif_id) |>
        summarize(across(everything(), ~ paste0(unique(.x), collapse = ";")), .groups = "drop"),
      by = "motif_id"
    ) |>
    relocate(tf_family, .after = motif_id) |>
    select(-any_of("motif_alt_id")) |>
    relocate(c(tf_id, tf_name), .after = tf_family)
}

collapse_unique_chr <- function(x, sep = ";") {
  values <- x |>
    as.character() |>
    str_trim() |>
    discard(~ is.na(.x) || .x == "")

  if (length(values) == 0) {
    return(NA_character_)
  }

  paste(sort(unique(values)), collapse = sep)
}

build_fimo_ref_alt_comparison <- function(fimo_hits_at_variant, variant_tbl) {
  if (nrow(fimo_hits_at_variant) == 0) {
    return(variant_tbl[0, ] |>
      transmute(
        variant_id,
        rsid,
        motif_id = character(),
        tf_family = character(),
        tf_id = character(),
        tf_name = character(),
        database = character(),
        db = character(),
        start = integer(),
        end = integer(),
        width = integer(),
        strand = character(),
        matched_sequence_ref = character(),
        matched_sequence_alt = character(),
        score_ref = numeric(),
        score_alt = numeric(),
        pvalue_ref = numeric(),
        pvalue_alt = numeric(),
        qvalue_ref = numeric(),
        qvalue_alt = numeric(),
        delta_llr = numeric(),
        delta_llr_abs = numeric(),
        best_allele = character(),
        best_score = numeric(),
        best_pvalue = numeric(),
        best_qvalue = numeric(),
        has_positive_score = logical(),
        passes_significance_filter = logical(),
        passes_priority_filter = logical()
      ))
  }

  fimo_hits_at_variant |>
    mutate(
      variant_id = coalesce(variant_id, str_remove(seqnames, "_(ref|alt)$")),
      allele = coalesce(allele, str_extract(seqnames, "(ref|alt)$"))
    ) |>
    group_by(variant_id, start, end, width, strand, motif_id, tf_family, tf_id, tf_name, database, db) |>
    filter(n() == 2, setequal(allele, c("ref", "alt"))) |>
    summarize(
      matched_sequence_ref = matched_sequence[allele == "ref"][1],
      matched_sequence_alt = matched_sequence[allele == "alt"][1],
      score_ref = score[allele == "ref"][1],
      score_alt = score[allele == "alt"][1],
      pvalue_ref = pvalue[allele == "ref"][1],
      pvalue_alt = pvalue[allele == "alt"][1],
      qvalue_ref = qvalue[allele == "ref"][1],
      qvalue_alt = qvalue[allele == "alt"][1],
      .groups = "drop"
    ) |>
    mutate(
      delta_llr = score_alt - score_ref,
      delta_llr_abs = abs(delta_llr),
      best_allele = if_else(score_ref >= score_alt, "ref", "alt"),
      best_score = pmax(score_ref, score_alt, na.rm = TRUE),
      best_pvalue = if_else(best_allele == "ref", pvalue_ref, pvalue_alt),
      best_qvalue = if_else(best_allele == "ref", qvalue_ref, qvalue_alt),
      has_positive_score = score_ref > 0 | score_alt > 0,
      passes_significance_filter = (!is.na(best_qvalue) & best_qvalue < 0.05) |
        (!is.na(best_pvalue) & best_pvalue < 1e-04),
      passes_priority_filter = has_positive_score & passes_significance_filter
    ) |>
    left_join(variant_tbl |> select(variant_id, rsid), by = "variant_id") |>
    relocate(rsid, .after = variant_id) |>
    arrange(desc(delta_llr_abs), best_pvalue, best_qvalue)
}

build_tf_support_summary <- function(interval_hits) {
  if (nrow(interval_hits) == 0) {
    return(tibble(
      variant_id = character(),
      tf_name = character(),
      supported_in_remap = logical(),
      supported_in_chip_atlas = logical(),
      supported_in_dnase_hint16 = logical(),
      supported_in_dnase_hint20 = logical(),
      any_interval_support = logical(),
      matched_sources = character(),
      chip_atlas_cell_types = character(),
      chip_atlas_tissue_types = character(),
      chip_atlas_feature_labels = character(),
      remap_cell_types = character(),
      remap_feature_labels = character(),
      dnase_hint16_motifs = character(),
      dnase_hint20_motifs = character()
    ))
  }

  interval_hits |>
    transmute(
      variant_id,
      source_dataset,
      cell_type = as.character(cell_type),
      tissue_type = as.character(tissue_type),
      feature_label = as.character(feature_label),
      gene = na_if(as.character(gene), ""),
      motif = na_if(as.character(motif), "")
    ) |>
    pivot_longer(
      cols = c(gene, motif),
      names_to = "match_field",
      values_to = "tf_name"
    ) |>
    filter(!is.na(tf_name), tf_name != "") |>
    distinct() |>
    group_by(variant_id, tf_name) |>
    summarize(
      supported_in_remap = any(source_dataset == "ReMap"),
      supported_in_chip_atlas = any(source_dataset == "ChIP-Atlas"),
      supported_in_dnase_hint16 = any(source_dataset == "DNase_HINT16"),
      supported_in_dnase_hint20 = any(source_dataset == "DNase_HINT20"),
      any_interval_support = n() > 0,
      matched_sources = collapse_unique_chr(source_dataset),
      chip_atlas_cell_types = collapse_unique_chr(cell_type[source_dataset == "ChIP-Atlas"]),
      chip_atlas_tissue_types = collapse_unique_chr(tissue_type[source_dataset == "ChIP-Atlas"]),
      chip_atlas_feature_labels = collapse_unique_chr(feature_label[source_dataset == "ChIP-Atlas"]),
      remap_cell_types = collapse_unique_chr(cell_type[source_dataset == "ReMap"]),
      remap_feature_labels = collapse_unique_chr(feature_label[source_dataset == "ReMap"]),
      dnase_hint16_motifs = collapse_unique_chr(feature_label[source_dataset == "DNase_HINT16"]),
      dnase_hint20_motifs = collapse_unique_chr(feature_label[source_dataset == "DNase_HINT20"]),
      .groups = "drop"
    )
}

run_fimo <- function(rsSeqs, motifs, meme_path, thresh = 1, bfile = NULL, rds_path = NULL, force = FALSE) {
  if (!is.null(rds_path) && file.exists(rds_path) && !force) {
    cached_res <- readRDS(rds_path)
    cached_seqnames <- tryCatch(as.character(cached_res$seqnames), error = function(e) character())
    expected_seqnames <- names(rsSeqs)

    if (length(cached_seqnames) > 0 && all(cached_seqnames %in% expected_seqnames)) {
      log_message("Loading cached FIMO result from ", rds_path)
      return(cached_res)
    }

    log_message(
      "Cached FIMO result at ", rds_path,
      " does not match the current sequence IDs. Re-running FIMO and overwriting the cache."
    )
  }

  log_message("Running FIMO on ", length(rsSeqs), " sequences using ", basename(motifs))

  res <- runFimo(
    rsSeqs,
    text = FALSE,
    motifs = motifs,
    meme_path = meme_path,
    thresh = thresh,
    bfile = bfile
  )

  if (!is.null(rds_path)) {
    saveRDS(res, rds_path)
  }

  res
}

find_interval_overlaps_tbl <- function(query_tbl, target_tbl, target_chr_col, target_start_col, target_end_col, target_dataset_name) {
  empty_result <- bind_cols(
    query_tbl[0, ],
    target_tbl[0, ] |>
      mutate(
        dataset_chr = character(),
        source_dataset = character()
      )
  )

  target_filtered <- target_tbl |>
    mutate(
      dataset_chr = normalize_chr(.data[[target_chr_col]])
    ) |>
    semi_join(query_tbl |> distinct(chr), by = c("dataset_chr" = "chr"))

  if (nrow(target_filtered) == 0) {
    return(empty_result)
  }

  query_gr <- make_snp_granges(query_tbl)
  target_gr <- GRanges(
    seqnames = format_chr_for_ucsc(target_filtered$dataset_chr),
    ranges = IRanges(
      start = as.numeric(target_filtered[[target_start_col]]) + 1L * (target_start_col == "chromStart"),
      end = as.numeric(target_filtered[[target_end_col]])
    )
  )

  hits <- findOverlaps(query_gr, target_gr)

  if (length(hits) == 0) {
    return(empty_result)
  }

  bind_cols(
    query_tbl[queryHits(hits), ],
    target_filtered[subjectHits(hits), ] |> mutate(source_dataset = target_dataset_name)
  )
}

merged_regulatory_features_columns <- c(
  "chrom",
  "start",
  "end",
  "source_dataset",
  "source_record_id",
  "experiment_type",
  "cell_type",
  "tissue_type",
  "gene",
  "motif",
  "feature_label",
  "biosample_group",
  "score",
  "strand",
  "metadata_compact"
)

empty_merged_regulatory_features_overlap_tbl <- function(query_tbl) {
  out <- query_tbl[0, ]

  character_cols <- c(
    "chrom",
    "source_dataset",
    "source_record_id",
    "experiment_type",
    "cell_type",
    "tissue_type",
    "gene",
    "motif",
    "feature_label",
    "biosample_group",
    "strand",
    "metadata_compact"
  )

  numeric_cols <- c("start", "end", "score")

  for (col in character_cols) {
    out[[col]] <- character()
  }

  for (col in numeric_cols) {
    out[[col]] <- numeric()
  }

  out
}

read_merged_regulatory_features_hits_for_chr <- function(query_chr_tbl, bgz_dir, file_prefix = "merged_regulatory_features") {
  chr_value <- unique(query_chr_tbl$chr)

  if (length(chr_value) != 1) {
    stop("read_merged_regulatory_features_hits_for_chr expects a single chromosome at a time.")
  }

  bgz_path <- file.path(bgz_dir, paste0(file_prefix, "_chr", chr_value, ".tsv.bgz"))
  tbi_path <- paste0(bgz_path, ".tbi")

  if (!file.exists(bgz_path)) {
    log_message("Merged regulatory features file not found for chr", chr_value, ": ", bgz_path)
    return(tibble())
  }

  if (!file.exists(tbi_path)) {
    stop("Missing tabix index for merged regulatory features file: ", tbi_path)
  }

  query_positions <- sort(unique(query_chr_tbl$pos))
  query_gr <- GRanges(
    seqnames = format_chr_for_ucsc(chr_value),
    ranges = IRanges(start = query_positions, end = query_positions)
  )

  tbx <- TabixFile(bgz_path, index = tbi_path)
  on.exit(close(tbx), add = TRUE)
  open(tbx)

  lines_by_region <- scanTabix(tbx, param = query_gr)

  hit_tbls <- map2(
    lines_by_region,
    query_positions,
    function(lines, query_pos) {
      if (length(lines) == 0) {
        return(NULL)
      }

      query_chr <- chr_value

      dt <- data.table::fread(
        text = paste(lines, collapse = "\n"),
        sep = "\t",
        header = FALSE,
        quote = "",
        col.names = merged_regulatory_features_columns,
        na.strings = c("", "NA"),
        showProgress = FALSE
      )

      dt[, start := suppressWarnings(as.integer(start))]
      dt[, end := suppressWarnings(as.integer(end))]
      dt[, score := suppressWarnings(as.numeric(score))]
      dt[, query_chr := query_chr]
      dt[, query_pos := query_pos]
      as_tibble(dt)
    }
  ) |>
    compact()

  if (length(hit_tbls) == 0) {
    return(tibble())
  }

  bind_rows(hit_tbls)
}

find_merged_regulatory_features_overlaps_tbi <- function(query_tbl, bgz_dir, file_prefix = "merged_regulatory_features") {
  hits_tbl <- query_tbl |>
    distinct(chr, pos) |>
    group_split(chr, .keep = TRUE) |>
    map_dfr(
      ~ read_merged_regulatory_features_hits_for_chr(
        query_chr_tbl = .x,
        bgz_dir = bgz_dir,
        file_prefix = file_prefix
      )
    )

  if (nrow(hits_tbl) == 0) {
    return(empty_merged_regulatory_features_overlap_tbl(query_tbl))
  }

  query_tbl |>
    inner_join(
      hits_tbl,
      by = c("chr" = "query_chr", "pos" = "query_pos")
    )
}

write_summary_yaml <- function(config, path) {
  lines <- c(
    paste0("run_id: ", config$run_id),
    paste0("output_dir: ", config$output_dir),
    paste0("variant_count: ", config$variant_count),
    paste0("chromosomes: ", paste(config$chromosomes, collapse = ",")),
    paste0("run_fimo: ", tolower(as.character(config$run_fimo))),
    paste0("flank_bp: ", config$flank_bp),
    paste0("custom_dataset_enabled: ", tolower(as.character(config$custom_dataset_enabled))),
    paste0("custom_dataset_name: ", config$custom_dataset_name)
  )

  writeLines(lines, con = path)
}

infer_chr_from_path <- function(path) {
  filename <- basename(path)
  match <- str_match(filename, "(?i)(?:^|[^[:alnum:]])chr?([0-9]+|x|y|m|mt)(?:[^[:alnum:]]|$)")
  chr <- match[, 2]

  ifelse(is.na(chr), NA_character_, normalize_chr(chr))
}

resolve_custom_dataset_files <- function(single_file = NULL, file_list = NULL, file_glob = NULL) {
  files <- character()

  if (!is.null(single_file)) {
    files <- c(files, single_file)
  }

  if (!is.null(file_list)) {
    files <- c(
      files,
      str_split(file_list, ",", simplify = TRUE) |>
        as.character() |>
        str_trim()
    )
  }

  if (!is.null(file_glob)) {
    globbed <- Sys.glob(file_glob)
    files <- c(files, globbed)
  }

  files <- unique(files[nzchar(files)])

  if (length(files) == 0) {
    return(character())
  }

  missing <- files[!file.exists(files)]
  if (length(missing) > 0) {
    stop("Custom dataset file(s) not found: ", paste(missing, collapse = ", "))
  }

  files
}

read_tabular_dataset <- function(path, delim = "auto", has_header = TRUE) {
  chosen_delim <- delim
  if (chosen_delim == "auto") {
    chosen_delim <- if (str_detect(path, "\\.csv(\\.gz)?$")) "csv" else "tsv"
  }

  if (chosen_delim == "csv") {
    readr::read_csv(path, show_col_types = FALSE, col_names = has_header)
  } else {
    readr::read_tsv(path, show_col_types = FALSE, col_names = has_header)
  }
}

load_custom_dataset <- function(files, dataset_name, chr_col = NULL, start_col, end_col, delim = "auto", has_header = TRUE) {
  if (length(files) == 0) {
    return(NULL)
  }

  dataset_tbl <- map_dfr(
    files,
    function(path) {
      tbl <- read_tabular_dataset(path, delim = delim, has_header = has_header)
      tbl$source_file <- path
      tbl$file_chr_inferred <- infer_chr_from_path(path)
      tbl
    }
  )

  if (is.null(chr_col)) {
    if (all(!is.na(dataset_tbl$file_chr_inferred))) {
      dataset_tbl <- dataset_tbl |>
        mutate(custom_chr = file_chr_inferred)
      chr_col <- "custom_chr"
    } else {
      stop(
        "Custom dataset chromosome column was not provided and chromosome could not be inferred from all filenames. ",
        "Supply --custom-dataset-chr-col or use filenames containing values like chr1, chr2, chrX."
      )
    }
  }

  required_cols <- c(chr_col, start_col, end_col)
  missing_cols <- setdiff(required_cols, colnames(dataset_tbl))

  if (length(missing_cols) > 0) {
    stop(
      "Custom dataset is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  dataset_tbl |>
    mutate(
      custom_dataset_name = dataset_name,
      custom_chr_normalized = normalize_chr(.data[[chr_col]]),
      custom_start_numeric = as.numeric(.data[[start_col]]),
      custom_end_numeric = as.numeric(.data[[end_col]])
    )
}

###############################################################################
# Pipeline config
###############################################################################

output_dir <- get_arg("output-dir", file.path(getwd(), paste0("pipeline_output_", timestamp_string())))
run_id <- get_arg("run-id", basename(output_dir))
variants_file <- get_arg("variants-file", config_aliases = c("variants_file", "variant_file"))
variants_arg <- get_arg("variants")
run_fimo_flag <- coerce_flag(get_arg("run-fimo"), default = TRUE)
force_fimo_flag <- coerce_flag(get_arg("force-fimo"), default = FALSE)
flank_bp <- as.integer(get_arg("flank-bp", 20))
custom_dataset_name <- get_arg("custom-dataset-name", "custom_dataset")
custom_dataset_file <- get_arg("custom-dataset-file")
custom_dataset_files <- get_arg("custom-dataset-files")
custom_dataset_glob <- get_arg("custom-dataset-glob")
custom_dataset_chr_col <- get_arg("custom-dataset-chr-col")
custom_dataset_start_col <- get_arg("custom-dataset-start-col")
custom_dataset_end_col <- get_arg("custom-dataset-end-col")
custom_dataset_delim <- coerce_delim(get_arg("custom-dataset-delim"), default = "auto")
custom_dataset_has_header <- coerce_flag(get_arg("custom-dataset-has-header"), default = TRUE)
custom_dataset_zero_based <- coerce_flag(get_arg("custom-dataset-zero-based"), default = FALSE)

paths <- default_paths
for (path_name in names(default_paths)) {
  cli_name <- gsub("_", "-", path_name)
  config_name <- path_name
  paths[[path_name]] <- get_arg(cli_name, default_paths[[path_name]], config_aliases = c(config_name))
}
paths$merged_regulatory_features_bgz_dir <- get_arg(
  "merged-regulatory-features-bgz-dir",
  paths$merged_regulatory_features_bgz_dir,
  config_aliases = c("merged_regulatory_features_bgz_dir")
)
paths$merged_regulatory_features_prefix <- get_arg(
  "merged-regulatory-features-prefix",
  paths$merged_regulatory_features_prefix,
  config_aliases = c("merged_regulatory_features_prefix")
)

path_base_dirs <- unique(c(
  if (!is.null(config_path_cli)) dirname(config_path_cli) else character(),
  getwd(),
  project_root
))

for (path_name in c("jaspar_motifs", "jaspar_transfac", "cisbp_meme", "cisbp_info", "fimo_bg")) {
  if (!is.null(paths[[path_name]])) {
    paths[[path_name]] <- resolve_input_path(paths[[path_name]], base_dirs = path_base_dirs, path_type = "file")
  }
}

if (!is.null(paths$meme_bin)) {
  paths$meme_bin <- resolve_input_path(paths$meme_bin, base_dirs = path_base_dirs, path_type = "dir")
}

paths$merged_regulatory_features_bgz_dir <- resolve_input_path(
  paths$merged_regulatory_features_bgz_dir,
  base_dirs = path_base_dirs,
  path_type = "dir"
)

if (!is.null(variants_file)) {
  variants_file <- resolve_input_path(variants_file, base_dirs = path_base_dirs, path_type = "file")
}

if (!is.null(custom_dataset_file)) {
  custom_dataset_file <- resolve_input_path(custom_dataset_file, base_dirs = path_base_dirs, path_type = "file")
}

if (!is.null(custom_dataset_files)) {
  custom_dataset_files <- map_chr(
    str_split(custom_dataset_files, ",")[[1]],
    ~ resolve_input_path(str_trim(.x), base_dirs = path_base_dirs, path_type = "file")
  )
}

validate_pipeline_inputs <- function(paths, run_fimo_flag) {
  if (!dir.exists(paths$merged_regulatory_features_bgz_dir)) {
    stop(
      "Merged regulatory features bundle directory not found: ", paths$merged_regulatory_features_bgz_dir,
      ". Set merged_regulatory_features_bgz_dir in the config or with --merged-regulatory-features-bgz-dir."
    )
  }

  if (run_fimo_flag) {
    required_files <- c(
      jaspar_motifs = paths$jaspar_motifs,
      jaspar_transfac = paths$jaspar_transfac,
      cisbp_meme = paths$cisbp_meme,
      cisbp_info = paths$cisbp_info
    )

    missing_files <- names(required_files)[vapply(required_files, function(x) is.null(x) || !file.exists(x), logical(1))]

    if (length(missing_files) > 0) {
      stop(
        "run_fimo=true but required motif resource file(s) are missing: ",
        paste(missing_files, collapse = ", "),
        ". Provide them in the config or place them under the repo data/ directory."
      )
    }

    if (is.null(paths$meme_bin) || !dir.exists(paths$meme_bin) || !file.exists(file.path(paths$meme_bin, "fimo"))) {
      stop(
        "run_fimo=true but meme_bin does not contain the MEME Suite 'fimo' executable: ",
        ifelse(is.null(paths$meme_bin), "NULL", paths$meme_bin),
        ". Set meme_bin in the config or ensure 'fimo' is on PATH before running the pipeline."
      )
    }

    if (!is.null(paths$fimo_bg) && !file.exists(paths$fimo_bg)) {
      warning(
        "FIMO background model not found at ", paths$fimo_bg,
        ". Proceeding without a background model."
      )
      paths$fimo_bg <- NULL
    }
  }

  paths
}

paths <- validate_pipeline_inputs(paths, run_fimo_flag)

custom_dataset_input_files <- resolve_custom_dataset_files(
  single_file = custom_dataset_file,
  file_list = custom_dataset_files,
  file_glob = custom_dataset_glob
)

custom_dataset_enabled <- length(custom_dataset_input_files) > 0

if (custom_dataset_enabled && (is.null(custom_dataset_start_col) || is.null(custom_dataset_end_col))) {
  stop(
    "Custom dataset input requires --custom-dataset-start-col and --custom-dataset-end-col."
  )
}

pipeline_dirs <- list(
  output = ensure_dir(output_dir),
  results = ensure_dir(file.path(output_dir, "results")),
  tables = ensure_dir(file.path(output_dir, "results", "tables")),
  rds = ensure_dir(file.path(output_dir, "results", "rds")),
  logs = ensure_dir(file.path(output_dir, "logs")),
  diagnostics = ensure_dir(file.path(output_dir, "diagnostics")),
  config = ensure_dir(file.path(output_dir, "config"))
)

log_file <- file.path(pipeline_dirs$logs, "pipeline.log")
pipeline_log_con <- file(log_file, open = "at")

on.exit({
  if (!is.null(pipeline_log_con) && isOpen(pipeline_log_con)) {
    close(pipeline_log_con)
  }
}, add = TRUE)

log_message("Starting pipeline run: ", run_id)
log_message("Output directory: ", output_dir)
log_message("Built-in reference and interval resources are expected to use hg38 coordinates.")

###############################################################################
# Pipeline
###############################################################################

tryCatch({
variant_tbl <- load_variants(
  variants_arg = variants_arg,
  variants_file = variants_file
)

variant_tbl <- get_rsids(variant_tbl)

config_tbl <- list(
  run_id = run_id,
  output_dir = output_dir,
  variant_count = nrow(variant_tbl),
  chromosomes = sort(unique(variant_tbl$chr)),
  run_fimo = run_fimo_flag,
  flank_bp = flank_bp,
  custom_dataset_enabled = custom_dataset_enabled,
  custom_dataset_name = if (custom_dataset_enabled) custom_dataset_name else "none"
)

write_summary_yaml(config_tbl, file.path(pipeline_dirs$config, "run_config.yml"))
safe_write_tsv(variant_tbl, file.path(pipeline_dirs$tables, "input_variants.tsv.gz"))

diagnostics <- list()

diagnostics[["variant_summary"]] <- variant_tbl |>
  count(chr, name = "n_variants") |>
  arrange(chr)

###############################################################################
# Motif metadata
###############################################################################

motif_meta <- NULL

if (run_fimo_flag) {
  log_message("Building motif metadata tables")

  motif_meta <- build_motif_info(
    cisbp_info_filepath = paths$cisbp_info,
    jaspar_transfac_filepath = paths$jaspar_transfac
  )

  safe_write_tsv(motif_meta$motif_info, file.path(pipeline_dirs$tables, "motif_metadata.tsv.gz"))
  safe_save_rds(motif_meta, file.path(pipeline_dirs$rds, "motif_metadata.rds"))
} else {
  log_message("Skipping motif metadata build because run_fimo=false")
}

###############################################################################
# Merged regulatory features bundle overlaps
###############################################################################

log_message("Querying merged tabix-indexed regulatory features bundle")

merged_regulatory_features_hits <- find_merged_regulatory_features_overlaps_tbi(
  query_tbl = variant_tbl,
  bgz_dir = paths$merged_regulatory_features_bgz_dir,
  file_prefix = paths$merged_regulatory_features_prefix
)

remap_hits <- merged_regulatory_features_hits |>
  filter(source_dataset == "ReMap") |>
  mutate(
    name = source_record_id,
    tf_name = gene
  )

chip_atlas_hits <- merged_regulatory_features_hits |>
  filter(source_dataset == "ChIP-Atlas") |>
  filter(experiment_type != "histone_chipseq")

dnase_hint16_hits <- merged_regulatory_features_hits |>
  filter(source_dataset == "DNase_HINT16")

dnase_hint20_hits <- merged_regulatory_features_hits |>
  filter(source_dataset == "DNase_HINT20")

safe_write_tsv(remap_hits, file.path(pipeline_dirs$tables, "remap_overlaps.tsv.gz"))
safe_write_tsv(chip_atlas_hits, file.path(pipeline_dirs$tables, "chip_atlas_overlaps.tsv.gz"))
safe_write_tsv(dnase_hint16_hits, file.path(pipeline_dirs$tables, "dnase_hint16_overlaps.tsv.gz"))
safe_write_tsv(dnase_hint20_hits, file.path(pipeline_dirs$tables, "dnase_hint20_overlaps.tsv.gz"))

diagnostics[["remap_summary"]] <- remap_hits |>
  count(chr, name = "n_overlaps") |>
  arrange(chr)

diagnostics[["chip_atlas_summary"]] <- chip_atlas_hits |>
  count(chr, name = "n_overlaps") |>
  arrange(chr)

diagnostics[["dnase_hint16_summary"]] <- dnase_hint16_hits |>
  count(chr, name = "n_overlaps") |>
  arrange(chr)

diagnostics[["dnase_hint20_summary"]] <- dnase_hint20_hits |>
  count(chr, name = "n_overlaps") |>
  arrange(chr)

###############################################################################
# Custom user-supplied dataset overlaps
###############################################################################

custom_dataset_hits <- tibble()

if (custom_dataset_enabled) {
  log_message("Reading custom dataset input: ", custom_dataset_name)

  custom_dataset_tbl <- load_custom_dataset(
    files = custom_dataset_input_files,
    dataset_name = custom_dataset_name,
    chr_col = custom_dataset_chr_col,
    start_col = custom_dataset_start_col,
    end_col = custom_dataset_end_col,
    delim = custom_dataset_delim,
    has_header = custom_dataset_has_header
  )

  custom_dataset_start_overlap_col <- "custom_start_numeric"
  custom_dataset_end_overlap_col <- "custom_end_numeric"

  if (custom_dataset_zero_based) {
    custom_dataset_tbl <- custom_dataset_tbl |>
      mutate(custom_start_numeric = custom_start_numeric + 1)
  }

  custom_dataset_hits <- find_interval_overlaps_tbl(
    query_tbl = variant_tbl,
    target_tbl = custom_dataset_tbl,
    target_chr_col = "custom_chr_normalized",
    target_start_col = custom_dataset_start_overlap_col,
    target_end_col = custom_dataset_end_overlap_col,
    target_dataset_name = custom_dataset_name
  )

  safe_write_tsv(
    custom_dataset_hits,
    file.path(pipeline_dirs$tables, paste0(custom_dataset_name, "_overlaps.tsv.gz"))
  )

  diagnostics[[paste0(custom_dataset_name, "_summary")]] <- custom_dataset_hits |>
    count(chr, name = "n_overlaps") |>
    arrange(chr)
}

###############################################################################
# Sequence construction and FIMO
###############################################################################

if (run_fimo_flag) {
  allele_seq_tbl <- build_allele_sequences(variant_tbl, flank = flank_bp)

  variant_levels <- unique(allele_seq_tbl$variant_id)
  allele_seq_tbl <- allele_seq_tbl |>
    mutate(
      variant_index = match(variant_id, variant_levels),
      seq_id = paste0("variant", variant_index, "_", allele)
    )

  safe_write_tsv(allele_seq_tbl, file.path(pipeline_dirs$tables, "allele_sequences.tsv.gz"))

  sequence_diagnostics <- allele_seq_tbl |>
    distinct(variant_id, center_base_original, ref) |>
    mutate(matches_reference = center_base_original == ref)

  diagnostics[["sequence_reference_check"]] <- sequence_diagnostics

  log_message("Preparing allele-specific sequences for FIMO")

  rsSeqs_all <- DNAStringSet(allele_seq_tbl$updated_sequence)
  names(rsSeqs_all) <- allele_seq_tbl$seq_id

  seq_metadata <- allele_seq_tbl |>
    select(seqnames = seq_id, variant_id, rsid, chr, pos, allele, ref, alt)

  fimo_jaspar_rds <- file.path(pipeline_dirs$rds, "fimo_jaspar_results.rds")
  fimo_cisbp_rds <- file.path(pipeline_dirs$rds, "fimo_cisbp_results.rds")

  fimo_results_jaspar <- run_fimo(
    rsSeqs = rsSeqs_all,
    motifs = paths$jaspar_motifs,
    meme_path = paths$meme_bin,
    thresh = 1,
    bfile = paths$fimo_bg,
    rds_path = fimo_jaspar_rds,
    force = force_fimo_flag
  )

  fimo_results_cisbp <- run_fimo(
    rsSeqs = rsSeqs_all,
    motifs = paths$cisbp_meme,
    meme_path = paths$meme_bin,
    thresh = 1,
    bfile = paths$fimo_bg,
    rds_path = fimo_cisbp_rds,
    force = force_fimo_flag
  )

  fimo_results_all <- combine_fimo_results(
    fimo_results_cisbp = fimo_results_cisbp,
    fimo_results_jaspar = fimo_results_jaspar,
    rsSeqs = rsSeqs_all,
    motif_info_long = motif_meta$motif_info_long,
    seq_metadata = seq_metadata
  )

  fimo_hits_at_variant <- fimo_results_all |>
    filter(start <= flank_bp + 1, end >= flank_bp + 1)

  fimo_pairs <- fimo_hits_at_variant |>
    mutate(
      variant_id = coalesce(variant_id, str_remove(seqnames, "_(ref|alt)$")),
      allele = coalesce(allele, str_extract(seqnames, "(ref|alt)$"))
    ) |>
    group_by(variant_id, start, end, width, strand, motif_id, tf_family, tf_id, tf_name, database, db) |>
    filter(n() == 2, setequal(allele, c("ref", "alt"))) |>
    ungroup()

  tf_support_summary <- build_tf_support_summary(merged_regulatory_features_hits)

  fimo_ref_alt_comparison <- build_fimo_ref_alt_comparison(
    fimo_hits_at_variant = fimo_hits_at_variant,
    variant_tbl = variant_tbl
  ) |>
    left_join(tf_support_summary, by = c("variant_id", "tf_name")) |>
    mutate(
      across(
        c(
          supported_in_remap,
          supported_in_chip_atlas,
          supported_in_dnase_hint16,
          supported_in_dnase_hint20,
          any_interval_support
        ),
        ~ replace_na(.x, FALSE)
      )
    )

  fimo_ref_alt_prioritized <- fimo_ref_alt_comparison |>
    filter(passes_priority_filter) |>
    arrange(
      desc(any_interval_support),
      desc(supported_in_chip_atlas),
      desc(delta_llr_abs),
      best_pvalue,
      best_qvalue
    )

  safe_write_tsv(fimo_hits_at_variant, file.path(pipeline_dirs$tables, "fimo_hits_at_variant.tsv.gz"))
  safe_write_tsv(fimo_pairs, file.path(pipeline_dirs$tables, "fimo_hits_ref_alt_pairs.tsv.gz"))
  safe_write_tsv(
    fimo_ref_alt_comparison,
    file.path(pipeline_dirs$tables, "fimo_ref_alt_comparison.tsv.gz")
  )
  safe_write_tsv(
    fimo_ref_alt_prioritized,
    file.path(pipeline_dirs$results, "fimo_ref_alt_prioritized.tsv.gz")
  )
  safe_save_rds(fimo_results_all, file.path(pipeline_dirs$rds, "fimo_combined_results.rds"))

  diagnostics[["fimo_summary"]] <- tibble(
    metric = c(
      "n_fimo_hits_at_variant",
      "n_ref_alt_paired_hits",
      "n_variants_with_paired_hits",
      "n_ref_alt_comparisons",
      "n_prioritized_ref_alt_comparisons"
    ),
    value = c(
      nrow(fimo_hits_at_variant),
      nrow(fimo_pairs),
      n_distinct(fimo_pairs$variant_id),
      nrow(fimo_ref_alt_comparison),
      nrow(fimo_ref_alt_prioritized)
    )
  )
}

###############################################################################
# Diagnostics and session info
###############################################################################

diagnostics_combined <- imap_dfr(
  diagnostics,
  ~ mutate(.x, diagnostic_name = .y, .before = 1)
)

safe_write_tsv(diagnostics_combined, file.path(pipeline_dirs$diagnostics, "pipeline_diagnostics.tsv.gz"))

summary_tbl <- tibble(
  metric = c(
    "n_input_variants",
    "n_input_chromosomes",
    "n_remap_hits",
    "n_chip_atlas_hits",
    "n_dnase_hint16_hits",
    "n_dnase_hint20_hits",
    "n_custom_dataset_hits"
  ),
  value = c(
    nrow(variant_tbl),
    n_distinct(variant_tbl$chr),
    nrow(remap_hits),
    nrow(chip_atlas_hits),
    nrow(dnase_hint16_hits),
    nrow(dnase_hint20_hits),
    nrow(custom_dataset_hits)
  )
)

safe_write_tsv(summary_tbl, file.path(pipeline_dirs$results, "run_summary.tsv"))

session_info_path <- file.path(pipeline_dirs$diagnostics, "sessionInfo.txt")
capture.output(sessionInfo(), file = session_info_path)

log_message("Pipeline complete")
}, error = function(e) {
  error_text <- conditionMessage(e)
  log_message("Pipeline failed: ", error_text)
  message("Pipeline failed: ", error_text)
  message("See log: ", log_file)
  quit(save = "no", status = 1)
})
