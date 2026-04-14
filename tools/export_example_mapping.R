args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript tools/export_example_mapping.R <input_rda> <output_csv>")
}

input_rda <- args[[1]]
output_csv <- args[[2]]

env <- new.env(parent = emptyenv())
load(input_rda, envir = env)
if (!exists("peptidomics_protein_mapping_example", envir = env)) {
  stop("The RDA file does not contain peptidomics_protein_mapping_example")
}

mapping <- get("peptidomics_protein_mapping_example", envir = env)
write.csv(mapping, output_csv, row.names = FALSE)

