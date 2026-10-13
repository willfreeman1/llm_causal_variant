# Run the published CALDERA model on the 1,348 UK Biobank neighborhoods.
# One call per trait, because a gene's score depends on the trait and CALDERA
# matches scores by gene id only.

args <- commandArgs(trailingOnly = TRUE)
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/"))
root <- normalizePath(file.path(script_dir, ".."), winslash = "/")
trait <- NULL
if (length(args) >= 1 && dir.exists(args[[1]])) {
  root <- normalizePath(args[[1]], winslash = "/")
  if (length(args) >= 2) trait <- args[[2]]
} else if (length(args) >= 1) {
  trait <- args[[1]]
}
caldera_path <- file.path(root, "third_party/caldera_run")
indir <- file.path(root, "data/ukb/caldera_inputs")
outdir <- file.path(root, "data/ukb/caldera_results")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

source(file.path(caldera_path, "z_caldera.R"))

cs_paths <- list.files(indir, pattern = "_cs\\.tsv$", full.names = TRUE)
traits <- sub("_cs\\.tsv$", "", basename(cs_paths))
if (!is.null(trait)) {
  keep <- traits == trait
  if (!any(keep)) stop("Unknown trait: ", trait)
  cs_paths <- cs_paths[keep]
  traits <- traits[keep]
}

out <- vector("list", length(traits))
t0 <- proc.time()
for (i in seq_along(traits)) {
  message(date(), " trait ", i, "/", length(traits), ": ", traits[i])
  genes <- caldera(
    pops_file = file.path(indir, paste0(traits[i], "_pops.tsv")),
    cs_file = cs_paths[i],
    assembly = 37,
    caldera_path = caldera_path
  )
  genes$trait <- traits[i]
  out[[i]] <- genes
}
elapsed <- proc.time() - t0
message("elapsed sec ", elapsed[["elapsed"]])

result <- data.table::rbindlist(out)
outfile <- file.path(outdir, "weeks_caldera_genes.tsv")
data.table::fwrite(result, outfile, sep = "\t")
message("wrote ", outfile, " rows ", nrow(result), " loci ", length(unique(result$locus)))
