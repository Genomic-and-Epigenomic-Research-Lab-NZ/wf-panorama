#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("Usage: methylcibersort.R <input_beta_matrix> <mixture_matrix_name> <base_sig_matrix_filename> <sample_name> [cancer_type]")

input_beta_matrix <- args[1]
mixture_matrix_name <- args[2]
base_sig_matrix_filename <- args[3]
sample_name <- args[4]
cancer_type <- args[5]

# helper to read CSV/TSV
read_input <- function(path) {
    tryCatch(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
        error = function(e) {
            tryCatch(read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE),
                error = function(e2) stop(paste("Failed to read input beta matrix:", path))
            )
        }
    )
}

cat("Reading beta matrix from:", input_beta_matrix, "\n")
beta_df <- read_input(input_beta_matrix)

# # Attempt to coerce to probe x sample matrix
# if ("probe" %in% names(beta_df) && "beta" %in% names(beta_df)) {
#     mix_df <- beta_df[, c("probe", "beta")]
#     names(mix_df) <- c("NAME", sample_name)
# } else if ("probe" %in% names(beta_df) && ncol(beta_df) > 1) {
#     mix_df <- beta_df
#     names(mix_df)[1] <- "NAME"
#     # if only two columns, rename second to sample
#     if (ncol(mix_df) == 2) names(mix_df)[2] <- sample_name
# } else if (ncol(beta_df) >= 2) {
#     names(beta_df)[1] <- "NAME"
#     mix_df <- beta_df
#     if (ncol(mix_df) == 2) names(mix_df)[2] <- sample_name
# } else {
#     stop("Unrecognised beta matrix format")
# }

rownames(beta_df) <- beta_df$probe
beta_df$probe <- NULL
colnames(beta_df) <- sample_name
beta_df <- as.matrix(beta_df)

# # write mixture matrix as tab-delimited
# mix_out <- paste0(mixture_matrix_name, ".txt")
# cat("Writing mixture matrix to:", mix_out, "\n")
# write.table(
#     mix_df,
#     file = mix_out, sep = "\t", row.names = FALSE, quote = FALSE
# )

# If MethylCIBERSORT is available, attempt to use its helper; otherwise write signature file placeholder
# suppressWarnings(suppressMessages({
#     have_mc <- requireNamespace("MethylCIBERSORT", quietly = TRUE)
# }))

library(MethylCIBERSORT)

cat("MethylCIBERSORT available; preparing signature and mixture using package functions\n")
# load signatures
data("V2_Signatures")
sig_key <- paste0(cancer_type, "_v2_Signature.txt")
cat("Using cancer type signature:", sig_key, "\n")
if (!sig_key %in% names(Signatures)) {
    stop(paste(
        "Cancer type signature not found in V2_Signatures:", sig_key,
        "\nAvailable signatures:", paste(names(Signatures), collapse = ", ")
    ))
}
base_sig_matrix <- Signatures[[sig_key]]

# Export signature to file
cat("Writing minimal signature to: ", base_sig_matrix_filename, "\n")
write.table(
    base_sig_matrix,
    file = base_sig_matrix_filename,
    sep = "\t", row.names = FALSE, quote = FALSE
)

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
Prep.CancerType(
    Beta = beta_df,
    Probes = base_sig_matrix$NAME,
    fname = mixture_matrix_name
)

# # If MethylCIBERSORT's Prep.CancerType is available, call it (best-effort)
# if (exists("Prep.CancerType", where = asNamespace("MethylCIBERSORT"), inherits = FALSE)) {
#     tryCatch(
#         {
#             MethylCIBERSORT::Prep.CancerType(Beta = as.matrix(mix_df[, -1, drop = FALSE]), Probes = mix_df$NAME, fname = mixture_matrix_name)
#             cat("Called Prep.CancerType to prepare mixture files\n")
#         },
#         error = function(e) {
#             cat("Prep.CancerType failed:", e$message, "\n")
#         }
#     )
# }


cat("methylcibersort.R completed.\n")
