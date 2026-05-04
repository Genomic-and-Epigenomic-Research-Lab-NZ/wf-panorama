FROM ubuntu:22.04

# Copy files into the image
COPY mCS_venv_test.yaml /tmp/mCS_venv_test.yaml
COPY MethylCIBERSORT_0.2.0.tar.gz /opt/MethylCIBERSORT_0.2.0.tar.gz

# Install system dependencies
RUN apt-get update && apt-get install -y \
        wget \
        curl \
        ca-certificates \
        build-essential \
        gfortran \
        git \
        bzip2 \
        libarchive-dev \
    && rm -rf /var/lib/apt/lists/*

# Install Mambaforge
RUN curl -L https://github.com/conda-forge/mambaforge/releases/latest/download/Mambaforge-Linux-x86_64.sh -o /tmp/mambaforge.sh \
    && bash /tmp/mambaforge.sh -b -p /opt/mambaforge \
    && rm /tmp/mambaforge.sh

# Set environment variables
ENV PATH="/opt/mambaforge/bin:$PATH"

# Default entrypoint
ENTRYPOINT []
CMD ["/bin/bash"]
