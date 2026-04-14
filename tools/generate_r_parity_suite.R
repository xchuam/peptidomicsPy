args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript tools/generate_r_parity_suite.R <upstream_dir> <fixtures_dir> <output_dir>")
}

upstream_dir <- normalizePath(args[[1]], mustWork = TRUE)
fixtures_dir <- normalizePath(args[[2]], mustWork = TRUE)
output_dir <- normalizePath(args[[3]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

lib_user <- Sys.getenv("R_LIBS_USER")
if (nzchar(lib_user)) {
  dir.create(lib_user, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(lib_user, .libPaths()))
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(limma)
  library(patchwork)
  library(gridExtra)
  library(cowplot)
  library(ggseqlogo)
})

theme_pubr <- function(base_size = 12, base_family = "",
                       border = FALSE, margin = TRUE,
                       legend = c("top", "bottom", "left", "right", "none"),
                       x.text.angle = 0) {
  half_line <- base_size / 2
  if (!is.numeric(legend)) legend <- match.arg(legend)
  if (x.text.angle > 5) xhjust <- 1 else xhjust <- NULL

  if (border) {
    panel.border <- element_rect(fill = NA, colour = "black", linewidth = 0.7)
    axis.line <- element_blank()
  } else {
    panel.border <- element_blank()
    axis.line <- element_line(colour = "black", linewidth = 0.5)
  }

  if (margin) {
    plot.margin <- margin(half_line, half_line, half_line, half_line)
  } else {
    plot.margin <- unit(c(0.5, 0.3, 0.3, 0.3), "mm")
  }

  .theme <- theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      panel.border = panel.border,
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = axis.line,
      axis.text = element_text(color = "black"),
      legend.key = element_blank(),
      strip.background = element_rect(fill = "#F2F2F2", colour = "black", linewidth = 0.7),
      plot.margin = plot.margin,
      legend.position = legend,
      complete = TRUE
    )

  if (x.text.angle != 0) {
    .theme <- .theme + theme(axis.text.x = element_text(angle = x.text.angle, hjust = xhjust))
  }

  .theme
}

write_gz_csv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
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

write_plot_png <- function(plot_obj, path, width, height, dpi = 150) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(filename = path, plot = plot_obj, width = width, height = height, dpi = dpi, limitsize = FALSE)
}

r_files <- list.files(file.path(upstream_dir, "R"), pattern = "\\.R$", full.names = TRUE)
invisible(lapply(r_files, source))

result <- processPeptides(
  peptides_file = file.path(fixtures_dir, "Yogurtexample_QR188-205.csv"),
  intensity_columns_file = file.path(fixtures_dir, "Intensity_columns.csv"),
  protein_mapping_file = file.path(fixtures_dir, "protein_mapping.csv")
)

dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "plots"), recursive = TRUE, showWarnings = FALSE)

# Table parity cases
filter_exact <- filterPeptides(
  result = result,
  seqs = c("AAGGPGAPADPGRPT", "AAIDEASKKLNAQ")
)
filter_regex <- filterPeptides(
  result = result,
  seq_pattern = "^AA"
)
filter_grouped <- filterPeptides(
  result = result,
  seq_pattern = "DQAM",
  filter_params = list(Yogurt = "Y1", Digest.stage = "G120")
)
filter_union <- filterPeptides(
  result = result,
  seqs = c("AAGGPGAPADPGRPT"),
  seq_pattern = "^AAI"
)

write_filter_case <- function(name, x) {
  out_dir <- file.path(output_dir, "tables", "filter", name)
  write_gz_csv(x$dt.peptides, file.path(out_dir, "dt_peptides.csv.gz"))
  write_gz_csv(x$dt.peptides.int, file.path(out_dir, "dt_peptides_int.csv.gz"))
  write_gz_csv(x$dt.peptides.int.reps, file.path(out_dir, "dt_peptides_int_reps.csv.gz"))
  write_gz_csv(x$dt.peptides.typenum, file.path(out_dir, "dt_peptides_typenum.csv.gz"))
  write_gz_csv(x$dt.peptides.typenum.reps, file.path(out_dir, "dt_peptides_typenum_reps.csv.gz"))
  write_gz_csv(x$dt.int_col, file.path(out_dir, "dt_int_col.csv.gz"))
}

write_filter_case("exact", filter_exact)
write_filter_case("regex", filter_regex)
write_filter_case("grouped", filter_grouped)
write_filter_case("union", filter_union)

ttest_plain_main <- ttestPeptides(
  result = result,
  comparisons = list(c("G120_Y1", "I120_Y1")),
  test_method = "plain"
)
ttest_plain_reordered <- ttestPeptides(
  result = result,
  comparisons = list(c("Y1_G120", "Y1_I120")),
  test_method = "plain"
)
ttest_plain_explicit <- ttestPeptides(
  result = result,
  comparisons = list(c("Digest.stage=G120_Yogurt=Y1", "Digest.stage=I120_Yogurt=Y1")),
  test_method = "plain"
)
ttest_plain_pooled <- ttestPeptides(
  result = result,
  comparisons = list(c("Y1", "Y2")),
  test_method = "plain"
)
ttest_plain_equalvar <- ttestPeptides(
  result = result,
  comparisons = list(c("G120_Y1", "I120_Y1")),
  test_method = "plain",
  equal_var = TRUE
)
ttest_treat_main <- ttestPeptides(
  result = result,
  comparisons = list(c("G120_Y1", "I120_Y1")),
  test_method = "treat"
)

write_gz_csv(ttest_plain_main[[1]], file.path(output_dir, "tables", "ttest", "plain_main.csv.gz"))
write_gz_csv(ttest_plain_reordered[[1]], file.path(output_dir, "tables", "ttest", "plain_reordered.csv.gz"))
write_gz_csv(ttest_plain_explicit[[1]], file.path(output_dir, "tables", "ttest", "plain_explicit.csv.gz"))
write_gz_csv(ttest_plain_pooled[[1]], file.path(output_dir, "tables", "ttest", "plain_pooled.csv.gz"))
write_gz_csv(ttest_plain_equalvar[[1]], file.path(output_dir, "tables", "ttest", "plain_equalvar.csv.gz"))
write_gz_csv(ttest_treat_main[[1]], file.path(output_dir, "tables", "ttest", "treat_main.csv.gz"))

# Plot parity cases
write_plot_png(
  plot_int(
    result,
    type = "mean",
    x_var = "Yogurt",
    filter_params = list(Digest.stage = "G120"),
    color_by = "Protein.name"
  ),
  file.path(output_dir, "plots", "plot_int_mean_named.png"),
  width = 8,
  height = 5
)

write_plot_png(
  plot_int(result, type = "reps"),
  file.path(output_dir, "plots", "plot_int_reps_default.png"),
  width = 8,
  height = 5
)

write_plot_png(
  plot_type_num(
    result,
    type = "reps",
    x_var = "Yogurt",
    facet_rows = "Digest.stage",
    color_by = "Protein.group"
  ),
  file.path(output_dir, "plots", "plot_type_num_reps.png"),
  width = 8,
  height = 5
)

write_plot_png(
  plot_count(
    result,
    type = "mean",
    x_var = "Yogurt",
    color_by = "none"
  ),
  file.path(output_dir, "plots", "plot_count_mean_none.png"),
  width = 8,
  height = 5
)

write_plot_png(
  plot_length_distribution(
    result,
    metric = "intensity",
    filter_params = list(Digest.stage = "G120"),
    facet_rows = "Yogurt"
  ),
  file.path(output_dir, "plots", "plot_length_bar.png"),
  width = 8,
  height = 5
)

write_plot_png(
  plot_length_distribution(
    result,
    type = "reps",
    metric = "type_num",
    plot_mode = "density",
    facet_rows = "Digest.stage"
  ),
  file.path(output_dir, "plots", "plot_length_density.png"),
  width = 8,
  height = 5
)

write_plot_png(
  plot_gravy_vs_intensity(result),
  file.path(output_dir, "plots", "plot_gravy_scatter.png"),
  width = 8,
  height = 5
)

write_plot_png(
  plot_gravy_vs_intensity(
    result,
    type = "reps",
    plot_mode = "density",
    facet_rows = "Digest.stage",
    filter_params = list(Yogurt = "Y1")
  ),
  file.path(output_dir, "plots", "plot_gravy_density.png"),
  width = 8,
  height = 5
)

write_plot_png(
  plot_volcano(
    ttest_plain_main,
    comparisons = list(c("G120_Y1", "I120_Y1")),
    test_method = "plain"
  )[[1]],
  file.path(output_dir, "plots", "plot_volcano_plain.png"),
  width = 8,
  height = 5
)

treat_labels <- head(ttest_treat_main[[1]]$Sequence, 2)
write_plot_png(
  plot_volcano(
    ttest_treat_main,
    comparisons = list(c("G120_Y1", "I120_Y1")),
    test_method = "treat",
    label_seqs = treat_labels,
    highlight_seqs = list("black" = treat_labels[1])
  )[[1]],
  file.path(output_dir, "plots", "plot_volcano_treat_annotated.png"),
  width = 8,
  height = 5
)

write_plot_png(
  plot_pep_align(
    result,
    protein_name = "P81265",
    protein_seq = paste(rep("A", 650), collapse = ""),
    auto_size = FALSE
  ),
  file.path(output_dir, "plots", "plot_pep_align_mean.png"),
  width = 12,
  height = 8
)

write_plot_png(
  plot_cleavage_site(result, terminal = "N"),
  file.path(output_dir, "plots", "plot_cleavage_N.png"),
  width = 8,
  height = 5
)

write_plot_png(
  plot_cleavage_site(
    result,
    terminal = "both",
    measure = "type_num",
    replicate_mode = "reps",
    filter_params = list(Yogurt = "Y1")
  ),
  file.path(output_dir, "plots", "plot_cleavage_both_filtered.png"),
  width = 12,
  height = 5
)
