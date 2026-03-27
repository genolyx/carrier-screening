FROM mambaorg/micromamba:1.5.8

USER root

# Install system dependencies
RUN apt-get update && apt-get install -y \
    curl \
    wget \
    git \
    procps \
    openjdk-17-jre-headless \
    && rm -rf /var/lib/apt/lists/*

# Install Nextflow
RUN curl -s https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin/ && \
    chmod +x /usr/local/bin/nextflow && \
    nextflow -version

# Create application directories
WORKDIR /app
RUN mkdir -p /data/fastq /data/analysis /data/output /data/log \
    && mkdir -p /var/log/supervisor

# Copy conda environment file
COPY environment.yml .

# Create conda environment
RUN micromamba install -y -n base -f environment.yml && \
    micromamba clean --all --yes

# Copy entrypoint script
COPY docker/entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

# Set environment variables
ENV PATH="/opt/conda/bin:${PATH}"
ENV BASE_DIR=/data

ENTRYPOINT ["/entrypoint.sh"]
CMD ["bash"]
