#!/usr/bin/env Rscript
# ==============================================================================
# gRodon2 Maximum Growth Rate Estimation Wrapper for MAGs
# ==============================================================================
# Supports both single-MAG and batch mode (e.g. 50 MAGs per single R session):
#   Rscript scripts/run_grodon2.R --batch_file <batch.tsv> --out <batch_growth.tsv> [--temp 20] [--mode partial]
#   Rscript scripts/run_grodon2.R --batch_dir <ffn_dir> --out <batch_growth.tsv> [--temp 20] [--mode partial]
#   Rscript scripts/run_grodon2.R --ffn <path.ffn> --mag <mag_id> --out <output.tsv> [--temp 20] [--mode partial]
# ==============================================================================

suppressPackageStartupMessages({
  library(methods)
})

# 1. Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) > 0 && idx + 1 <= length(args)) {
    return(args[idx + 1])
  }
  return(default)
}

ffn_path   <- get_arg("--ffn")
mag_id     <- get_arg("--mag", "unknown_mag")
batch_file <- get_arg("--batch_file")
batch_dir  <- get_arg("--batch_dir")
ribo_path  <- get_arg("--ribo_file")
out_path   <- get_arg("--out")
temp_val   <- as.numeric(get_arg("--temp", "20"))
mode_val   <- get_arg("--mode", "partial")

if (is.null(out_path) || (is.null(ffn_path) && is.null(batch_file) && is.null(batch_dir))) {
  cat("ERROR: Must provide --out and at least one of --ffn, --batch_file, or --batch_dir.\n", file = stderr())
  quit(status = 1)
}

# Ensure output directory exists
out_dir <- dirname(out_path)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

# 2. Base FASTA reader fallback
read_fasta_simple <- function(file_path) {
  lines <- readLines(file_path)
  headers <- grep("^>", lines)
  if (length(headers) == 0) return(character(0))
  
  starts <- headers + 1
  ends <- c(headers[-1] - 1, length(lines))
  
  seqs <- vapply(seq_along(headers), function(i) {
    if (starts[i] > ends[i]) return("")
    paste(lines[starts[i]:ends[i]], collapse = "")
  }, character(1))
  
  names(seqs) <- sub("^>([^\t ]+).*", "\\1", lines[headers])
  return(seqs)
}

# 3. Detect gRodon2 / gRodon package once
grod_pkg <- if (requireNamespace("gRodon2", quietly = TRUE)) {
  "gRodon2"
} else if (requireNamespace("gRodon", quietly = TRUE)) {
  "gRodon"
} else {
  NULL
}

# 4. Parse ribosomal annotations if provided
ribo_map <- list()
if (!is.null(ribo_path) && file.exists(ribo_path) && file.info(ribo_path)$size > 0) {
  tryCatch({
    r_lines <- readLines(ribo_path)
    r_lines <- r_lines[nzchar(trimws(r_lines))]
    if (length(r_lines) > 0) {
      header_parts <- strsplit(r_lines[1], "\t")[[1]]
      has_header <- any(grepl("^(MAG|Gene_ID|Target)", header_parts, ignore.case = TRUE))
      data_lines <- if (has_header) r_lines[-1] else r_lines
      for (l in data_lines) {
        parts <- strsplit(l, "\t")[[1]]
        if (length(parts) >= 2) {
          m <- trimws(parts[1])
          g <- trimws(parts[2])
          g_clean <- sub("^.*___", "", g)
          ribo_map[[m]] <- c(ribo_map[[m]], g, g_clean)
        } else if (length(parts) == 1) {
          g <- trimws(parts[1])
          if (grepl("___", g)) {
            m <- sub("___.*$", "", g)
            g_clean <- sub("^.*___", "", g)
            ribo_map[[m]] <- c(ribo_map[[m]], g, g_clean)
          }
        }
      }
      ribo_map <- lapply(ribo_map, unique)
      cat(sprintf("INFO: Loaded ribosomal annotations for %d MAGs from %s\n", length(ribo_map), ribo_path), file = stderr())
    }
  }, error = function(e) {
    cat(sprintf("WARNING: Failed parsing ribo_file %s: %s\n", ribo_path, e$message), file = stderr())
  })
}

# 5. Main prediction function for a single MAG
run_single_prediction <- function(ffn, mag, temp, mode_param, mag_ribos = NULL) {
  res_df <- data.frame(
    MAG = mag,
    DoublingTime_hours = NA_real_,
    MaxGrowthRate_per_hour = NA_real_,
    CUB_HEG = NA_real_,
    GrowthCategory = "Unknown",
    Method = ifelse(is.null(grod_pkg), "None", grod_pkg),
    Mode = mode_param,
    Temperature = temp,
    Status = "Pending",
    stringsAsFactors = FALSE
  )

  if (is.null(grod_pkg)) {
    res_df$Status <- "Package_Not_Installed (gRodon2/gRodon required in R environment)"
    return(res_df)
  }

  if (!file.exists(ffn) || file.info(ffn)$size == 0) {
    res_df$Status <- "Empty_or_Missing_FFN"
    return(res_df)
  }

  dna_seqs <- NULL
  if (requireNamespace("Biostrings", quietly = TRUE)) {
    tryCatch({
      dna_seqs <- Biostrings::readDNAStringSet(ffn)
    }, error = function(e) {
      dna_seqs <- NULL
    })
  }

  if (is.null(dna_seqs)) {
    simple_seqs <- read_fasta_simple(ffn)
    if (length(simple_seqs) == 0) {
      res_df$Status <- "No_Sequences_Found"
      return(res_df)
    }
  }

  tryCatch({
    if (is.null(dna_seqs) && requireNamespace("Biostrings", quietly = TRUE)) {
      dna_seqs <- Biostrings::readDNAStringSet(ffn)
    }

    seq_names <- names(dna_seqs)
    seq_ids <- sub("[ \t].*$", "", seq_names)

    # 1. Match from external ribosomal annotations (e.g. HMMER / CheckM)
    is_ribo_from_file <- rep(FALSE, length(seq_names))
    if (!is.null(mag_ribos) && length(mag_ribos) > 0) {
      is_ribo_from_file <- (seq_ids %in% mag_ribos) |
                           (seq_names %in% mag_ribos) |
                           (paste0(mag, "___", seq_ids) %in% mag_ribos)
    }

    # 2. Match from sequence header annotations (e.g. Bakta / Prokka / NCBI)
    is_ribo_from_header <- grepl("^(?!.*(methyl|hydroxy)).*0S ribosomal protein", seq_names, ignore.case = TRUE, perl = TRUE) |
                           grepl("ribosomal protein", seq_names, ignore.case = TRUE) |
                           grepl("\\b(rps|rpl|rpm)[A-Z0-9]+\\b", seq_names, ignore.case = TRUE)

    highly_expressed <- is_ribo_from_file | is_ribo_from_header
    num_ribos <- sum(highly_expressed)

    if (num_ribos < 10) {
      res_df$Status <- if (num_ribos == 0) {
        "No_Ribosomal_Proteins_Identified"
      } else {
        sprintf("Too_Few_Ribosomal_Proteins (found %d, min 10 required)", num_ribos)
      }
      return(res_df)
    }

    cat(sprintf("INFO [%s]: Running %s prediction (mode=%s, temp=%s°C, ribo_genes=%d)...\n", 
                mag, grod_pkg, mode_param, temp, num_ribos), file = stderr())

    pred_fn <- get("predictGrowth", asNamespace(grod_pkg))
    pred <- pred_fn(dna_seqs, highly_expressed = highly_expressed, mode = mode_param, temperature = temp)

    d_time <- as.numeric(pred$d)
    cub_val <- if (!is.null(pred$CUBHE)) as.numeric(pred$CUBHE) else (if (!is.null(pred$CUB)) as.numeric(pred$CUB) else NA_real_)

    res_df$DoublingTime_hours <- round(d_time, 3)
    res_df$MaxGrowthRate_per_hour <- round(log(2) / d_time, 4)
    res_df$CUB_HEG <- round(cub_val, 4)
    res_df$GrowthCategory <- ifelse(d_time < 5, "Copiotroph (Fast)", "Oligotroph (Slow)")
    res_df$Status <- "Success"
    return(res_df)
  }, error = function(e) {
    clean_err <- gsub("[\r\n\t\"]+", " ", as.character(e$message))
    clean_err <- trimws(clean_err)
    cat(sprintf("WARNING [%s]: %s prediction failed: %s\n", mag, grod_pkg, clean_err), file = stderr())
    res_df$Status <- paste0("Error: ", clean_err)
    return(res_df)
  })
}

# 6. Build list of jobs to process
jobs <- list()

if (!is.null(batch_file) && file.exists(batch_file)) {
  lines <- readLines(batch_file)
  lines <- lines[nzchar(trimws(lines))]
  for (line in lines) {
    parts <- strsplit(line, "[\t, ]+")[[1]]
    if (length(parts) >= 2) {
      jobs[[length(jobs) + 1]] <- list(mag = parts[1], ffn = parts[2])
    } else if (length(parts) == 1) {
      f_path <- parts[1]
      bname <- sub("\\.(ffn|fa|fasta|fna)(\\.gz)?$", "", basename(f_path))
      jobs[[length(jobs) + 1]] <- list(mag = bname, ffn = f_path)
    }
  }
} else if (!is.null(batch_dir) && dir.exists(batch_dir)) {
  ffn_files <- list.files(batch_dir, pattern = "\\.ffn(\\.gz)?$", full.names = TRUE, recursive = TRUE)
  for (f_path in ffn_files) {
    bname <- sub("\\.ffn(\\.gz)?$", "", basename(f_path))
    jobs[[length(jobs) + 1]] <- list(mag = bname, ffn = f_path)
  }
} else if (!is.null(ffn_path)) {
  jobs[[1]] <- list(mag = mag_id, ffn = ffn_path)
}

if (length(jobs) == 0) {
  cat("WARNING: No input FFN jobs found.\n", file = stderr())
  empty_df <- data.frame(
    MAG = character(0),
    DoublingTime_hours = numeric(0),
    MaxGrowthRate_per_hour = numeric(0),
    CUB_HEG = numeric(0),
    GrowthCategory = character(0),
    Method = character(0),
    Mode = character(0),
    Temperature = numeric(0),
    Status = character(0),
    stringsAsFactors = FALSE
  )
  write.table(empty_df, file = out_path, sep = "\t", quote = FALSE, row.names = FALSE)
  quit(status = 0)
}

# 7. Execute predictions sequentially in this single R process
results_list <- vector("list", length(jobs))
for (i in seq_along(jobs)) {
  job <- jobs[[i]]
  mag_ribos <- ribo_map[[job$mag]]
  results_list[[i]] <- run_single_prediction(job$ffn, job$mag, temp_val, mode_val, mag_ribos)
}

final_df <- do.call(rbind, results_list)

# Sanitize all string columns to guarantee clean TSV structure
for (col_name in names(final_df)) {
  if (is.character(final_df[[col_name]])) {
    final_df[[col_name]] <- gsub("[\r\n\t]+", " ", final_df[[col_name]])
    final_df[[col_name]] <- trimws(final_df[[col_name]])
  }
}

write.table(final_df, file = out_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("DONE: Wrote %d gRodon2 growth predictions to %s\n", nrow(final_df), out_path), file = stderr())
