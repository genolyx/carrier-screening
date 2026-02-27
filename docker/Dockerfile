FROM mambaorg/micromamba:1.5.8

USER root

# Install system dependencies
RUN apt-get update && apt-get install -y \
    curl \
    wget \
    git \
    procps \
    supervisor \
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

# Install Python packages for Dashboard and Daemon
RUN micromamba run -n base pip install --no-cache-dir \
    flask==3.0.0 \
    requests==2.31.0 \
    schedule==1.2.0 \
    watchdog==3.0.0 \
    gunicorn==21.2.0 \
    psutil

# Copy application files
COPY . /app/

# Copy supervisor configuration
COPY docker/supervisord.conf /etc/supervisor/conf.d/supervisord.conf

# Create entrypoint script
COPY docker/entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

# Expose ports
EXPOSE 5000 8080

# Set environment variables
ENV PATH="/opt/conda/bin:${PATH}"
ENV PYTHONUNBUFFERED=1
ENV BASE_DIR=/data

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=60s --retries=3 \
    CMD curl -f http://localhost:8080/api/health && curl -f http://localhost:5000/health || exit 1

# Run supervisor
ENTRYPOINT ["/entrypoint.sh"]
CMD ["/usr/bin/supervisord", "-c", "/etc/supervisor/supervisord.conf"]
