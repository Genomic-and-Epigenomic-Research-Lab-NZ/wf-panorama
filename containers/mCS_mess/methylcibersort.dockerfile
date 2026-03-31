FROM ubuntu:20.04

LABEL maintainer="Lucy Picard"
LABEL version="container to run R scripts with version 3 and MethylCIBERSORT"
LABEL description="R 3.6.3 + BioC + MethylCIBERSORT"

ENV R_LIBS_USER=/usr/local/lib/R/site-library
ENV PATH=/usr/local/bin:$PATH
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        wget curl git \
        gfortran \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        libblas-dev \
        liblapack-dev \
        libjpeg-dev \
        libtiff5-dev \
        libpng-dev \
        zlib1g-dev \
        libfreetype6-dev \
        libfontconfig1-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libbz2-dev \
        liblzma-dev \
        xz-utils \
        cmake \
        pkg-config \
        ca-certificates \
        python3 \
        python3-pip \
        libreadline-dev \
        libncurses5-dev && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Build and install R 3.6.3 from source
RUN set -eux; \
    cd /tmp; \
    wget -q https://cran.r-project.org/src/base/R-3/R-3.6.3.tar.gz -O R-3.6.3.tar.gz; \
    tar xzf R-3.6.3.tar.gz; \
    cd R-3.6.3; \
    ./configure --enable-R-shlib --with-blas --with-lapack --with-x=no --prefix=/usr/local; \
    make -j"$(nproc)"; \
    make install; \
    cd /; \
    rm -rf /tmp/R-3.6.3*; \
    R --version

RUN mkdir -p ${R_LIBS_USER}

# Add required resource files into the image
# These paths are relative to the build context; ensure the files exist at these locations when building.
COPY MethylCIBERSORT_0.2.0.tar.gz /opt/MethylCIBERSORT_0.2.0.tar.gz
COPY old_snkmk_mCS_env.txt /opt/old_snkmk_mCS_env.txt

# Create and run an R script to install CRAN/Bioconductor packages with versions from the env file
COPY install_R_packages.R /tmp/install_R_packages.R
RUN set -eux; \
    echo "Running R installer script to install packages (this may take a while)"; \
    Rscript --vanilla --verbose /tmp/install_R_packages.R

# Diagnostic: list installed R packages and save to csv for debugging
RUN set -eux; \
    Rscript -e "ip <- installed.packages(); write.csv(ip[,c('Package','Version')], '/tmp/installed_R_pkgs.csv', row.names=FALSE); print(head(ip[,c('Package','Version')],100))"; \
    cat /tmp/installed_R_pkgs.csv || true

# Install MethylCIBERSORT tarball
RUN R CMD INSTALL /opt/MethylCIBERSORT_0.2.0.tar.gz

# Default workdir
WORKDIR /home

# Minimal entrypoint
CMD ["R"]
