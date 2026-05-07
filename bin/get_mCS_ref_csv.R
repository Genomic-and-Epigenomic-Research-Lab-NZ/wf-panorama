#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

cancer_type <- args[1]
output_path <- args[2]

library(MethylCIBERSORT)
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
cat("Writing immune reference signature to: ", output_path, "\n")
write.table(
    base_sig_matrix,
    file = output_path,
    sep = ",", row.names = FALSE, quote = FALSE
)
