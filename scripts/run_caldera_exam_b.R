# CALDERA on Exam B credible sets. Assembly 38. Does not compare picks to labels.

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/"))
args_root <- normalizePath(file.path(script_dir, ".."), winslash = "/")
caldera_path <- file.path(args_root, "third_party/caldera_run")
indir <- file.path(args_root, "data/ukb/final_exam/caldera_inputs")
outdir <- file.path(args_root, "data/ukb/final_exam")

source(file.path(args_root, "third_party/caldera_run/z_caldera.R"))

cs_paths <- list.files(indir, pattern = "_cs\\.tsv$", full.names = TRUE)
traits <- sub("_cs\\.tsv$", "", basename(cs_paths))
out <- vector("list", length(traits))
for (i in seq_along(traits)) {
  message(date(), " trait ", i, "/", length(traits), ": ", traits[i])
  genes <- caldera(
    pops_file = file.path(indir, paste0(traits[i], "_pops.tsv")),
    cs_file = cs_paths[i],
    assembly = 38,
    caldera_path = caldera_path
  )
  genes$trait <- traits[i]
  out[[i]] <- genes
}
result <- data.table::rbindlist(out)
outfile <- file.path(outdir, "caldera_genes.tsv")
data.table::fwrite(result, outfile, sep = "\t")
message("wrote ", outfile, " rows ", nrow(result), " loci ", length(unique(result$locus)))
