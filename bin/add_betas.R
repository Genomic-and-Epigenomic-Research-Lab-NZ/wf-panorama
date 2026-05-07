#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: add_betas.R <input_pre_beta> <output_post_beta>")
}

pre_beta <- args[1]
post_beta <- args[2]

suppressMessages(library(minfi))
suppressMessages(library(data.table))

# Read input (try CSV then tab)
read_input <- function(path) {
    tryCatch(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
            error = function(e) {
                tryCatch(
                    read.table(
                        path,
                        header = TRUE, 
                        sep = "\t", 
                        stringsAsFactors = FALSE, 
                        check.names = FALSE
                    ),
                    error = function(e2) stop(paste(
                        "Failed to read input pre_beta:", 
                        path))
                )
            })
}

cat('Reading pre-beta from:', pre_beta, "\n")
df <- read_input(pre_beta)

# Ensure required columns present
if (!all(c('N','X') %in% names(df))) {
    stop('Input pre_beta must contain columns named N and X')
}

# Determine probe/position column to use as rownames
probe_col <- NULL
if ("position" %in% names(df)) {
    probe_col <- "position"
} else if ("probe" %in% names(df)) {
    probe_col <- "probe"
} else if ("pos" %in% names(df)) {
    probe_col <- "pos"
}

if (is.null(probe_col)) {
    # if no explicit probe column, create rownames from row numbers
    rownames_vec <- seq_len(nrow(df))
} else {
    rownames_vec <- as.character(df[[probe_col]])
}

# compute unmeth (N - X)
df$unmeth <- as.numeric(df$N) - as.numeric(df$X)

# Build Meth and Unmeth matrices
Meth <- matrix(as.numeric(df$X), ncol = 1)
Unmeth <- matrix(as.numeric(df$unmeth), ncol = 1)

rownames(Meth) <- rownames(Unmeth) <- rownames_vec
colnames(Meth) <- colnames(Unmeth) <- 'Sample_1'

# Create MethylSet and compute beta values
mset <- MethylSet(Meth = Meth, Unmeth = Unmeth)
beta_values <- getBeta(mset)

# add beta column to dataframe
df$beta <- beta_values

# bullshit code that breaks shit
# # Preserve original ordering
# if (!is.null(probe_col)) {
#     df$beta <- as.numeric(beta_values[as.character(df[[probe_col]]), 1])
# } else {
#     df$beta <- as.numeric(beta_values[,1])
# }

# drop helper column
df$unmeth <- NULL

# write output
cat('Writing post-beta to:', post_beta, "\n")
write.csv(df, file = post_beta, quote = FALSE, row.names = FALSE)
cat('Done.\n')
