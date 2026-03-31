logfile <- '/tmp/install_R_packages.log'
cat('Starting package install log\n', file = logfile)

env_file <- '/opt/old_snkmk_mCS_env.txt'
cat('Reading env file:', env_file, '\n', file = logfile, append = TRUE)
lines <- tryCatch(readLines(env_file), error = function(e) { cat('Failed to read env file: ', e, '\n', file = logfile); stop(e) })

# Parse env file into pinned map (name -> version) and mark Bioconductor vs CRAN
tok <- grep('^(r-|bioconductor-)', lines, value = TRUE)
cat('Found', length(tok), 'matching lines\n', file = logfile, append = TRUE)
# Normalize whitespace
tok <- gsub('\t', ' ', tok)
tok <- gsub(' +', ' ', tok)
parts <- strsplit(tok, ' ')

pinned <- list()   # named list: pinned[[pkg]] = list(version=..., source='CRAN'|'BioC')
bioc_pkgs <- character()

for (p in parts) {
  if (length(p) >= 2) {
    label <- p[[1]]
    ver <- p[[2]]
    if (startsWith(label, 'bioconductor-')) {
      name <- sub('^bioconductor-', '', label)
      name <- gsub('-', '', name)
      pinned[[name]] <- list(version = ver, source = 'BioC')
      bioc_pkgs <- unique(c(bioc_pkgs, name))
      cat('Pinned Bioconductor:', name, ver, '\n', file = logfile, append = TRUE)
    } else if (startsWith(label, 'r-')) {
      name <- sub('^r-', '', label)
      # Try common name variants
      candidates <- unique(c(name, gsub('-', '.', name), gsub('-', '', name)))
      # choose the most likely name: prefer dotted version if available
      chosen <- candidates[1]
      pinned[[chosen]] <- list(version = ver, source = 'CRAN')
      cat('Pinned CRAN:', chosen, ver, 'candidates:', paste(candidates, collapse=','), '\n', file = logfile, append = TRUE)
    }
  }
}

# Setup repos and helper packages
options(repos = c(CRAN = 'https://cran.rstudio.com'))

ensure_pkg <- function(pk) {
  if (!requireNamespace(pk, quietly = TRUE)) install.packages(pk, repos = 'https://cran.rstudio.com')
}

ensure_pkg('remotes')
ensure_pkg('BiocManager')
ensure_pkg('tools')

# Set Bioconductor release for R 3.6.x
cat('Setting Bioconductor version to 3.10\n', file = logfile, append = TRUE)
tryCatch({ BiocManager::install(version = '3.10', ask = FALSE, update = FALSE) }, error = function(e) { cat('BiocManager version selection failed: ', e, '\n', file = logfile, append = TRUE); stop(e) })

# Install Bioconductor packages as a batch using BiocManager (this installs Bioc-side dependencies)
if (length(bioc_pkgs) > 0) {
  cat('Installing Bioconductor packages via BiocManager:', paste(bioc_pkgs, collapse=','), '\n', file = logfile, append = TRUE)
  tryCatch({
    BiocManager::install(bioc_pkgs, ask = FALSE, update = FALSE, dependencies = TRUE)
  }, error = function(e) { cat('Bioconductor bulk install failed: ', conditionMessage(e), '\n', file = logfile, append = TRUE); stop(e) })
}

# Helper to read DESCRIPTION fields from a local tarball without installing
read_description_from_tar <- function(tarball) {
  td <- tempdir()
  d <- utils::untar(tarball, list = TRUE)
  # find the top-level folder inside tar
  topdir <- strsplit(d[1], '/')[[1]][1]
  unpack_dir <- file.path(td, paste0('pkg_unpack_', as.integer(Sys.time())))
  dir.create(unpack_dir, recursive = TRUE)
  utils::untar(tarball, exdir = unpack_dir)
  desc_path <- file.path(unpack_dir, topdir, 'DESCRIPTION')
  if (!file.exists(desc_path)) return(list())
  desc <- read.dcf(desc_path)
  res <- list()
  if ('Depends' %in% colnames(desc)) res$Depends <- desc[1, 'Depends']
  if ('Imports' %in% colnames(desc)) res$Imports <- desc[1, 'Imports']
  if ('LinkingTo' %in% colnames(desc)) res$LinkingTo <- desc[1, 'LinkingTo']
  # cleanup
  unlink(unpack_dir, recursive = TRUE)
  return(res)
}

# Parse dependency fields (Depends/Imports/LinkingTo) into package names (no version constraints)
parse_deps <- function(field) {
  if (is.null(field) || is.na(field) || nchar(field) == 0) return(character())
  # remove R version constraints and split
  parts <- unlist(strsplit(field, ','))
  parts <- gsub('\n', ' ', parts)
  parts <- trimws(parts)
  # remove entries like 'R (>= 3.5.0)'
  parts <- parts[!grepl('^R\b', parts)]
  # strip version constraints in parentheses
  parts <- gsub('\s*\(.*?\)', '', parts)
  parts <- trimws(parts)
  parts[parts != '']
}

# Recursive installer for a CRAN package pinned to a specific version
install_cran_pinned <- function(pkg, version, installed_cache = list()) {
  key <- paste0(pkg, '_', version)
  if (key %in% installed_cache) return(invisible(TRUE))
  # if package already installed and version matches, skip
  if (pkg %in% rownames(installed.packages())) {
    inst_ver <- as.character(packageVersion(pkg))
    if (inst_ver == version) {
      cat('Package', pkg, version, 'already installed\n', file = logfile, append = TRUE)
      installed_cache <<- c(installed_cache, key)
      return(invisible(TRUE))
    }
  }
  # ensure pinned entry exists
  if (is.null(pinned[[pkg]]) || pinned[[pkg]]$source != 'CRAN') {
    stop(paste('Package', pkg, 'not pinned as CRAN package in env file; cannot safely install dependency'))
  }
  cat('Downloading source for', pkg, version, '\n', file = logfile, append = TRUE)
  tarball <- tryCatch({ remotes::download_version(pkg, version = version, repos = 'https://cran.rstudio.com') }, error = function(e) { cat('download_version failed for', pkg, version, ':', conditionMessage(e), '\n', file = logfile, append = TRUE); stop(e) })
  # read DESCRIPTION to discover dependencies
  desc <- read_description_from_tar(tarball)
  deps <- unique(c(parse_deps(desc$Depends), parse_deps(desc$Imports), parse_deps(desc$LinkingTo)))
  # install dependencies first using pinned versions
  if (length(deps) > 0) {
    for (d in deps) {
      # skip base packages
      if (d %in% rownames(installed.packages(priority = 'base'))) next
      if (!is.null(pinned[[d]])) {
        if (pinned[[d]]$source == 'CRAN') {
          install_cran_pinned(d, pinned[[d]]$version, installed_cache)
        } else if (pinned[[d]]$source == 'BioC') {
          # Bioconductor dependency: attempt to install via BiocManager (BiocManager was already configured to release 3.10)
          cat('Installing Bioconductor dependency', d, 'via BiocManager\n', file = logfile, append = TRUE)
          tryCatch({ BiocManager::install(d, ask = FALSE, update = FALSE, dependencies = TRUE) }, error = function(e) { cat('Failed Bioc dependency', d, conditionMessage(e), '\n', file = logfile, append = TRUE); stop(e) })
        }
      } else {
        stop(paste('Dependency', d, 'of', pkg, 'is not pinned in env file; refusing to install unpinned dependency'))
      }
    }
  }
  # now install the package itself from the tarball without letting remotes resolve dependencies
  cat('Installing CRAN package', pkg, version, 'from tarball\n', file = logfile, append = TRUE)
  tryCatch({
    remotes::install_version(pkg, version = version, repos = 'https://cran.rstudio.com', upgrade = 'never', dependencies = FALSE)
  }, error = function(e) { cat('Failed to install', pkg, version, ':', conditionMessage(e), '\n', file = logfile, append = TRUE); stop(e) })
  installed_cache <<- c(installed_cache, key)
  invisible(TRUE)
}

# Install all pinned CRAN packages using topological resolution via the recursive installer
cran_to_install <- names(Filter(function(x) x$source == 'CRAN', pinned))
cat('CRAN packages to install (count):', length(cran_to_install), '\n', file = logfile, append = TRUE)
for (pkg in cran_to_install) {
  ver <- pinned[[pkg]]$version
  cat('Processing pinned CRAN package', pkg, ver, '\n', file = logfile, append = TRUE)
  install_cran_pinned(pkg, ver)
}

# Final sanity check: list installed packages
ip <- installed.packages()
write.csv(ip[,c('Package','Version')], '/tmp/installed_R_pkgs.csv', row.names = FALSE)
cat('Installed packages written to /tmp/installed_R_pkgs.csv\n', file = logfile, append = TRUE)
cat('Finished package installation\n', file = logfile, append = TRUE)
