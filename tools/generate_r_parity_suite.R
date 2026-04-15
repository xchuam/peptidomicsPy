args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript tools/generate_r_parity_suite.R <upstream_dir> <fixtures_dir> <output_dir>")
}

upstream_dir <- normalizePath(args[[1]], mustWork = TRUE)
fixtures_dir <- normalizePath(args[[2]], mustWork = TRUE)
output_dir <- normalizePath(args[[3]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
parity_scope <- Sys.getenv("PEPTIDOMICSPY_PARITY_SCOPE", unset = "full")
focus_plots_only <- identical(parity_scope, "plot_focus")

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

write_plain_csv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(df, path, row.names = FALSE)
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

mapping_env <- new.env(parent = emptyenv())
load(file.path(upstream_dir, "data", "peptidomics_protein_mapping_example.rda"), envir = mapping_env)
mapping_example <- data.table::as.data.table(get("peptidomics_protein_mapping_example", envir = mapping_env))
if (!focus_plots_only) {
  write_gz_csv(mapping_example, file.path(output_dir, "tables", "helpers", "mapping_example_dataset.csv.gz"))
}

gravy_examples <- data.table::data.table(
  peptide = c("ACDEFGHIK", "ALWKTML", "GGGGGGG", "", "AX"),
  gravy = vapply(c("ACDEFGHIK", "ALWKTML", "GGGGGGG", "", "AX"), calculate_gravy, numeric(1))
)
if (!focus_plots_only) {
  write_gz_csv(gravy_examples, file.path(output_dir, "tables", "helpers", "calculate_gravy_examples.csv.gz"))
}

write_process_case <- function(name, x) {
  out_dir <- file.path(output_dir, "tables", "process", name)
  write_gz_csv(x$dt.peptides, file.path(out_dir, "dt_peptides.csv.gz"))
  write_gz_csv(x$dt.peptides.int, file.path(out_dir, "dt_peptides_int.csv.gz"))
  write_gz_csv(x$dt.peptides.int.reps, file.path(out_dir, "dt_peptides_int_reps.csv.gz"))
  write_gz_csv(x$dt.peptides.int.ttest, file.path(out_dir, "dt_peptides_int_ttest.csv.gz"))
  write_gz_csv(x$dt.peptides.typenum, file.path(out_dir, "dt_peptides_typenum.csv.gz"))
  write_gz_csv(x$dt.peptides.typenum.reps, file.path(out_dir, "dt_peptides_typenum_reps.csv.gz"))
  write_gz_csv(x$dt.int_col, file.path(out_dir, "dt_int_col.csv.gz"))
  write_plain_csv(
    data.frame(
      grp_cols = paste(x$grp_cols, collapse = "|"),
      peptides_select_col_basic = paste(x$peptides_select_col.basic, collapse = "|")
    ),
    file.path(out_dir, "metadata.csv")
  )
}

if (!focus_plots_only) {
  write_process_case("default", result)
}

# Table parity cases
if (!focus_plots_only) {
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
    write_gz_csv(x$dt.peptides.int.ttest, file.path(out_dir, "dt_peptides_int_ttest.csv.gz"))
    write_gz_csv(x$dt.peptides.typenum, file.path(out_dir, "dt_peptides_typenum.csv.gz"))
    write_gz_csv(x$dt.peptides.typenum.reps, file.path(out_dir, "dt_peptides_typenum_reps.csv.gz"))
    write_gz_csv(x$dt.int_col, file.path(out_dir, "dt_int_col.csv.gz"))
  }

  write_filter_case("exact", filter_exact)
  write_filter_case("regex", filter_regex)
  write_filter_case("grouped", filter_grouped)
  write_filter_case("union", filter_union)
}

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
ttest_plain_custom_adjust <- ttestPeptides(
  result = result,
  comparisons = list(c("G120_Y1", "I120_Y1")),
  test_method = "plain",
  pseudocount = 5,
  adjust = "bonferroni"
)
ttest_plain_min_reps <- ttestPeptides(
  result = result,
  comparisons = list(c("G120_Y1", "I120_Y1")),
  test_method = "plain",
  min_reps_per_side = 3
)
ttest_treat_main <- ttestPeptides(
  result = result,
  comparisons = list(c("G120_Y1", "I120_Y1")),
  test_method = "treat"
)
ttest_treat_custom_threshold <- ttestPeptides(
  result = result,
  comparisons = list(c("G120_Y1", "I120_Y1")),
  test_method = "treat",
  lfc_thresh = 0.5,
  alpha = 0.1
)

write_gz_csv(ttest_plain_main[[1]], file.path(output_dir, "tables", "ttest", "plain_main.csv.gz"))
if (!focus_plots_only) {
  write_gz_csv(ttest_plain_reordered[[1]], file.path(output_dir, "tables", "ttest", "plain_reordered.csv.gz"))
  write_gz_csv(ttest_plain_explicit[[1]], file.path(output_dir, "tables", "ttest", "plain_explicit.csv.gz"))
  write_gz_csv(ttest_plain_pooled[[1]], file.path(output_dir, "tables", "ttest", "plain_pooled.csv.gz"))
  write_gz_csv(ttest_plain_equalvar[[1]], file.path(output_dir, "tables", "ttest", "plain_equalvar.csv.gz"))
  write_gz_csv(ttest_plain_custom_adjust[[1]], file.path(output_dir, "tables", "ttest", "plain_custom_adjust.csv.gz"))
  write_gz_csv(ttest_plain_min_reps[[1]], file.path(output_dir, "tables", "ttest", "plain_min_reps.csv.gz"))
}
write_gz_csv(ttest_treat_main[[1]], file.path(output_dir, "tables", "ttest", "treat_main.csv.gz"))
if (!focus_plots_only) {
  write_gz_csv(ttest_treat_custom_threshold[[1]], file.path(output_dir, "tables", "ttest", "treat_custom_threshold.csv.gz"))
}

# Plot parity cases
clone_result <- function(x) {
  list(
    dt.peptides = copy(x$dt.peptides),
    dt.peptides.int = copy(x$dt.peptides.int),
    dt.peptides.int.reps = copy(x$dt.peptides.int.reps),
    dt.peptides.int.ttest = copy(x$dt.peptides.int.ttest),
    dt.peptides.typenum = copy(x$dt.peptides.typenum),
    dt.peptides.typenum.reps = copy(x$dt.peptides.typenum.reps),
    dt.int_col = copy(x$dt.int_col),
    grp_cols = x$grp_cols,
    peptides_select_col.basic = x$peptides_select_col.basic
  )
}

result_pep_align_custom <- clone_result(result)
setnames(
  result_pep_align_custom$dt.peptides.int.reps,
  c("Sequence", "Start.position", "End.position", "Intensity"),
  c("Seq.custom", "Start.custom", "End.custom", "Intensity.custom")
)

pep_align_metadata <- data.frame(
  case = character(),
  suggested_fig_width = numeric(),
  suggested_fig_height = numeric(),
  save_file_location = character(),
  stringsAsFactors = FALSE
)

capture_pep_align_metadata <- function(case_name, plot_obj) {
  pep_align_metadata <<- rbind(
    pep_align_metadata,
    data.frame(
      case = case_name,
      suggested_fig_width = if (!is.null(attr(plot_obj, "suggested_fig_width"))) attr(plot_obj, "suggested_fig_width") else NA_real_,
      suggested_fig_height = if (!is.null(attr(plot_obj, "suggested_fig_height"))) attr(plot_obj, "suggested_fig_height") else NA_real_,
      save_file_location = if (!is.null(attr(plot_obj, "save_file_location"))) attr(plot_obj, "save_file_location") else "",
      stringsAsFactors = FALSE
    )
  )
}

if (!focus_plots_only) {
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
    plot_int(
      result,
      type = "reps",
      x_var = "Digest.stage",
      color_by = "Protein.group",
      facet_rows = "Yogurt",
      facet_cols = "Replicate",
      scientific_10_y = FALSE
    ),
    file.path(output_dir, "plots", "plot_int_faceted_linear.png"),
    width = 10,
    height = 6
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
    plot_type_num(
      result,
      type = "reps",
      x_var = "Digest.stage",
      filter_params = list(Yogurt = "Y1"),
      facet_cols = "Replicate",
      scientific_10_y = TRUE
    ),
    file.path(output_dir, "plots", "plot_type_num_filtered_scientific.png"),
    width = 9,
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
    plot_count(
      result,
      type = "reps",
      x_var = "Yogurt",
      color_by = "Protein.group",
      filter_params = list(Digest.stage = "G120"),
      facet_rows = "Digest.stage",
      facet_cols = "Replicate",
      scientific_10_y = TRUE
    ),
    file.path(output_dir, "plots", "plot_count_reps_faceted.png"),
    width = 10,
    height = 6
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
    plot_length_distribution(
      result,
      type = "reps",
      metric = "intensity",
      color_by = "Protein.name",
      facet_rows = "Digest.stage",
      facet_cols = "Yogurt+Replicate",
      scientific_10_y = FALSE
    ),
    file.path(output_dir, "plots", "plot_length_bar_named_linear.png"),
    width = 10,
    height = 6
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
    plot_gravy_vs_intensity(
      result,
      color_by = "Protein.name",
      facet_cols = "Digest.stage",
      alpha_value = 0.4
    ),
    file.path(output_dir, "plots", "plot_gravy_scatter_named.png"),
    width = 10,
    height = 5
  )
}

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
  plot_volcano(
    ttest_plain_main,
    comparisons = 1,
    test_method = "plain",
    show_threshold = TRUE,
    lfc_thresh = 0.5,
    alpha = 0.1,
    fill_values = c(no = "#AAAAAA", yes = "#E41A1C"),
    point_size = 2.5,
    point_alpha = 0.6,
    label_seqs = head(ttest_plain_main[[1]]$Sequence, 2),
    label_size = 4,
    label_col = "navy",
    highlight_seqs = setNames(
      list(head(ttest_plain_main[[1]]$Sequence, 1), head(ttest_plain_main[[1]]$Sequence, 2)[2]),
      c("black", "#377EB8")
    ),
    highlight_size = 3,
    highlight_stroke = 2,
    y_log_scale = TRUE
  )[[1]],
  file.path(output_dir, "plots", "plot_volcano_custom_ylog.png"),
  width = 8,
  height = 5
)

lookup_case_dir <- file.path(output_dir, "lookup_data_case")
dir.create(file.path(lookup_case_dir, "data"), recursive = TRUE, showWarnings = FALSE)
file.copy(
  file.path(upstream_dir, "data", "peptidomics_protein_mapping_example.rda"),
  file.path(lookup_case_dir, "data", "peptidomics_protein_mapping_example.rda"),
  overwrite = TRUE
)
old_wd <- getwd()
setwd(lookup_case_dir)
lookup_plot <- plot_pep_align(
  result,
  protein_name = "P02662"
)
setwd(old_wd)
capture_pep_align_metadata("plot_pep_align_lookup_default", lookup_plot)
write_plot_png(
  lookup_plot,
  file.path(output_dir, "plots", "plot_pep_align_lookup_default.png"),
  width = 8,
  height = 5
)

pep_align_mean_plot <- plot_pep_align(
  result,
  protein_name = "P81265",
  protein_seq = paste(rep("A", 650), collapse = ""),
  auto_size = FALSE
)
write_plot_png(
  pep_align_mean_plot,
  file.path(output_dir, "plots", "plot_pep_align_mean.png"),
  width = 12,
  height = 8
)

pep_align_ranged_plot <- plot_pep_align(
  result,
  protein_name = "P02662",
  filter_params = list(Digest.stage = "G120"),
  x_interval = 5,
  x_range = c(80, 165),
  y_range = c(0, 12),
  label_seq = list(highlight_1 = c("DQAMEDIKQ", "DQAMEDIKQM"), highlight_2 = c("AMEDIKQM")),
  label_col = list(highlight_1 = "red", highlight_2 = "blue")
)
write_plot_png(
  pep_align_ranged_plot,
  file.path(output_dir, "plots", "plot_pep_align_ranged_labels.png"),
  width = 9,
  height = 6
)

pep_align_custom_cols_plot <- plot_pep_align(
  result_pep_align_custom,
  protein_name = "P02662",
  type = "reps",
  filter_params = list(Yogurt = "Y1"),
  seq_col = "Seq.custom",
  start_col = "Start.custom",
  end_col = "End.custom",
  intensity_col = "Intensity.custom"
)
capture_pep_align_metadata("plot_pep_align_reps_custom_cols", pep_align_custom_cols_plot)
write_plot_png(
  pep_align_custom_cols_plot,
  file.path(output_dir, "plots", "plot_pep_align_reps_custom_cols.png"),
  width = 9,
  height = 6
)

pep_align_saved_path <- file.path(output_dir, "plots", "plot_pep_align_saved.png")
pep_align_saved_plot <- plot_pep_align(
  result,
  protein_name = "P02662",
  filter_params = list(Digest.stage = "G120"),
  auto_size = TRUE,
  save_file_location = pep_align_saved_path,
  save_file_dpi = 200,
  width_per_aa = 0.12,
  height_per_row = 0.2,
  min_width = 5,
  min_height = 4
)
capture_pep_align_metadata("plot_pep_align_saved", pep_align_saved_plot)

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

write_plot_png(
  plot_cleavage_site(
    result,
    terminal = "C",
    filter_params = list(Yogurt = "Y1", Digest.stage = "G120"),
    scientific_10_y = FALSE,
    drop_constant_groups = FALSE
  ),
  file.path(output_dir, "plots", "plot_cleavage_C_keep_constant_linear.png"),
  width = 8,
  height = 5
)

write_plain_csv(pep_align_metadata, file.path(output_dir, "tables", "helpers", "plot_pep_align_metadata.csv"))
