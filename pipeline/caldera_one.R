# One live CALDERA call. Writes a scratch file and does not replace the locked output.
args <- commandArgs(trailingOnly = TRUE)
root <- args[[1]]
cs_file <- args[[2]]
pops_file <- args[[3]]
out_file <- args[[4]]
caldera_path <- file.path(root, "third_party/caldera_run")
source(file.path(root, "third_party/caldera_run/z_caldera.R"))
genes <- caldera(
  pops_file = pops_file,
  cs_file = cs_file,
  assembly = 38,
  caldera_path = caldera_path
)
data.table::fwrite(genes, out_file, sep = "\t")
message("wrote ", out_file, " rows ", nrow(genes))
