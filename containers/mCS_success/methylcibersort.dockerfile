FROM ubuntu:22.04

# Copy files into the image
COPY conda-lock.yml /tmp/conda-lock.yml
COPY conda-lock.yml /opt/conda-lock.yml
COPY MethylCIBERSORT_0.2.0.tar.gz /opt/MethylCIBERSORT_0.2.0.tar.gz

# Install system dependencies
RUN apt-get update && apt-get install -y \
        wget \
        ca-certificates \
        build-essential \
        gfortran \
        git \
    && rm -rf /var/lib/apt/lists/*

# Install micromamba
RUN wget -qO- https://micromamba.snakepit.net/api/micromamba/linux-64/latest | tar -xvj -C /tmp/ \
    && mkdir -p /opt/micromamba/bin \
    && mv /tmp/bin/micromamba /opt/micromamba/bin/

# Create environment from lock file
RUN /opt/micromamba/bin/micromamba create -y --name mCS --file /tmp/conda-lock.yml

# Activate environment and install R package
RUN /opt/micromamba/bin/micromamba run -n mCS R CMD INSTALL /opt/MethylCIBERSORT_0.2.0.tar.gz

# Cleanup
RUN rm -rf /tmp/conda-lock.yml

# Set environment variables
ENV PATH="/opt/micromamba/envs/mCS/bin:/opt/micromamba/bin:$PATH"
ENV CONDA_PREFIX="/opt/micromamba/envs/mCS"

# Default entrypoint
ENTRYPOINT []
CMD ["/bin/bash"]
