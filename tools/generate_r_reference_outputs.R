args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript tools/generate_r_reference_outputs.R <upstream_dir> <fixtures_dir> <output_dir>")
}

upstream_dir <- normalizePath(args[[1]], mustWork = TRUE)
fixtures_dir <- normalizePath(args[[2]], mustWork = TRUE)
output_dir <- normalizePath(args[[3]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
})

r_files <- list.files(file.path(upstream_dir, "R"), pattern = "\\.R$", full.names = TRUE)
invisible(lapply(r_files, source))

write_gz_csv <- function(df, path) {
  for (col in names(df)) {
    if (is.list(df[[col]])) {
      df[[col]] <- vapply(
        df[[col]],
        function(x) {
          if (length(x) == 0 || all(is.na(x))) {
            return("")
          }
          paste(x, collapse = "|")
        },
        character(1)
      )
    }
  }
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  write.csv(df, con, row.names = FALSE)
}

result <- processPeptides(
  peptides_file = file.path(fixtures_dir, "Yogurtexample_QR188-205.csv"),
  intensity_columns_file = file.path(fixtures_dir, "Intensity_columns.csv"),
  protein_mapping_file = file.path(fixtures_dir, "protein_mapping.csv")
)

dir.create(file.path(output_dir, "processPeptides"), recursive = TRUE, showWarnings = FALSE)
write_gz_csv(result$dt.peptides, file.path(output_dir, "processPeptides", "dt_peptides.csv.gz"))
write_gz_csv(result$dt.peptides.int, file.path(output_dir, "processPeptides", "dt_peptides_int.csv.gz"))
write_gz_csv(result$dt.peptides.int.reps, file.path(output_dir, "processPeptides", "dt_peptides_int_reps.csv.gz"))
write_gz_csv(result$dt.peptides.int.ttest, file.path(output_dir, "processPeptides", "dt_peptides_int_ttest.csv.gz"))
write_gz_csv(result$dt.peptides.typenum, file.path(output_dir, "processPeptides", "dt_peptides_typenum.csv.gz"))
write_gz_csv(result$dt.peptides.typenum.reps, file.path(output_dir, "processPeptides", "dt_peptides_typenum_reps.csv.gz"))
write_gz_csv(result$dt.int_col, file.path(output_dir, "processPeptides", "dt_int_col.csv.gz"))

exact_filter <- filterPeptides(
  result = result,
  seqs = c("AAGGPGAPADPGRPT", "AAIDEASKKLNAQ")
)
regex_filter <- filterPeptides(
  result = result,
  seq_pattern = "^AA"
)
group_filter <- filterPeptides(
  result = result,
  seq_pattern = "DQAM",
  filter_params = list(Yogurt = "Y1", Digest.stage = "G120")
)

write_filter_outputs <- function(name, filter_result) {
  dir.create(file.path(output_dir, "filterPeptides", name), recursive = TRUE, showWarnings = FALSE)
  write_gz_csv(filter_result$dt.peptides, file.path(output_dir, "filterPeptides", name, "dt_peptides.csv.gz"))
  write_gz_csv(filter_result$dt.peptides.int, file.path(output_dir, "filterPeptides", name, "dt_peptides_int.csv.gz"))
  write_gz_csv(filter_result$dt.peptides.int.reps, file.path(output_dir, "filterPeptides", name, "dt_peptides_int_reps.csv.gz"))
  write_gz_csv(filter_result$dt.peptides.typenum, file.path(output_dir, "filterPeptides", name, "dt_peptides_typenum.csv.gz"))
  write_gz_csv(filter_result$dt.peptides.typenum.reps, file.path(output_dir, "filterPeptides", name, "dt_peptides_typenum_reps.csv.gz"))
  write_gz_csv(filter_result$dt.int_col, file.path(output_dir, "filterPeptides", name, "dt_int_col.csv.gz"))
}

write_filter_outputs("exact", exact_filter)
write_filter_outputs("regex", regex_filter)
write_filter_outputs("grouped", group_filter)

dir.create(file.path(output_dir, "ttestPeptides"), recursive = TRUE, showWarnings = FALSE)
ttest_plain <- ttestPeptides(
  result = result,
  comparisons = list(c("G120_Y1", "I120_Y1")),
  test_method = "plain"
)
ttest_treat <- ttestPeptides(
  result = result,
  comparisons = list(c("G120_Y1", "I120_Y1")),
  test_method = "treat"
)
write_gz_csv(ttest_plain[[1]], file.path(output_dir, "ttestPeptides", "plain_G120_Y1_vs_I120_Y1.csv.gz"))
write_gz_csv(ttest_treat[[1]], file.path(output_dir, "ttestPeptides", "treat_G120_Y1_vs_I120_Y1.csv.gz"))
