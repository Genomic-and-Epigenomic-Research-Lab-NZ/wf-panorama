#!/usr/bin/env Rscript

library(data.table)
library(MethylCIBERSORT)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("Usage: methylcibersort.R <input_beta_matrix> <mixture_matrix_name> <base_sig_matrix_filename> <sample_name> [cancer_type]")

input_beta_matrix <- args[1]
mixture_matrix_name <- args[2]
base_sig_matrix_filename <- args[3]
sample_name <- args[4]
cancer_type <- args[5]

# data("StromalMatrix_V2")
# table(Stromal_v2.pheno)

# input_beta_matrix <- "Test2.methatlas.csv"
# helper to read CSV/TSV - use for methatlas code
# read_input <- function(path) {
#     tryCatch(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
#         error = function(e) {
#             tryCatch(read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE),
#                 error = function(e2) stop(paste("Failed to read input beta matrix:", path))
#             )
#         }
#     )
# }

cat("Reading beta matrix from:", input_beta_matrix, "\n")
# input_beta_matrix <- "Test3_small.methatlas.csv"
beta_df <- fread(input_beta_matrix)

beta_df <- beta_df[, .(beta = mean(beta)), by = probe]

# Convert to data.frame
beta_df <- as.data.frame(beta_df)

# Make probes the index name
rownames(beta_df) <- beta_df$probe
# remove the probe column
beta_df$probe <- NULL
# retitle the values column with the sample
colnames(beta_df) <- sample_name
# convert to matrix for methylcibersort
beta_df <- as.matrix(beta_df)


cat("MethylCIBERSORT available; preparing signature and mixture using package functions\n")

# load signatures
data("V2_Signatures")
# cancer_type <- 'bladder'
sig_key <- paste0(cancer_type, "_v2_Signature.txt")
cat("Using cancer type signature:", sig_key, "\n")
if (!sig_key %in% names(Signatures)) {
    stop(paste(
        "Cancer type signature not found in V2_Signatures:", sig_key,
        "\nAvailable signatures:", paste(names(Signatures), collapse = ", ")
    ))
}
base_sig_matrix <- Signatures[[sig_key]]
dim(base_sig_matrix)
head(base_sig_matrix)

# Try the ref from methatlas
# methatlas_sig_matrix <- read_input('/projects/uow/GERL/dejlu879/Panorama/nanopore_multiBM_pipeline/resources/ref_atlas_bladder.csv')
# dim(methatlas_sig_matrix)
# head(methatlas_sig_matrix)

# adjust colnames to match the methylcibersort reference
# colnames(methatlas_sig_matrix) <- gsub("CpGs", "NAME", colnames(methatlas_sig_matrix))
# colnames(methatlas_sig_matrix) <- gsub("Bladder", "Cancer", colnames(methatlas_sig_matrix))
# colnames(methatlas_sig_matrix) <- gsub(sig_key, "Cancer", colnames(methatlas_sig_matrix))

# Export signature to file
cat("Writing minimal signature to: ", base_sig_matrix_filename, "\n")
write.table(
    base_sig_matrix,
    file = base_sig_matrix_filename,
    sep = "\t", row.names = FALSE, quote = FALSE
)
# base_sig_matrix_filename <- "Test2.bladder.mCS_ref.txt"
# write.table(
#     methatlas_sig_matrix,
#     file = base_sig_matrix_filename,
#     sep = "\t", row.names = FALSE, quote = FALSE
# )

# # Try to get signature object from package (best-effort)
# sig_out <- base_sig_matrix_file
# if (file.exists(sig_out)) {
#     cat("Base signature provided exists; copying to", sig_out, "\n")
#     file.copy(sig_out, sig_out, overwrite = TRUE)
# } else {
#     # write minimal signature (use probe names)
#     sig_df <- data.frame(NAME = mix_df$NAME, stringsAsFactors = FALSE)
#     write.table(sig_df, file = sig_out, sep = "\t", row.names = FALSE, quote = FALSE)
#     cat("Wrote minimal signature to:", sig_out, "\n")
# }

# below line was in the original code, not sure why
# load("/home/dejlu879/ProjectProtocol/bm_pipeline_dev/resources/mat.RData")

cat("Calling Prep.CancerType to prepare mixture files\n")
# mixture_matrix_name <- 'Test2.CS_mix_matrix'
    # Probes = methatlas_sig_matrix$NAME,
Prep.CancerType(
    Beta = beta_df,
    Probes = base_sig_matrix$NAME,
    fname = mixture_matrix_name
)

cat("methylcibersort.R completed.\n")
