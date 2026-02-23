FROM python:3.11-slim

LABEL maintainer="Dark Gene Pipeline Team"
LABEL description="Portal Integration Daemon for Dark Gene Pipeline"

# Install system dependencies
RUN apt-get update && apt-get install -y \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Create app directory
WORKDIR /app

# Copy requirements
COPY requirements.txt .

# Install Python dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Copy daemon files
COPY daemon.py .
COPY api_server.py .

# Create necessary directories
RUN mkdir -p /var/log /var/lib/daemon /pipeline

# Environment variables (override at runtime)
ENV PIPELINE_BASE_DIR=/pipeline
ENV PORTAL_URL=https://portal.genolyx.com
ENV PORTAL_API_KEY=""
ENV INSTITUTION_ID=default

# Expose API port
EXPOSE 8080

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD curl -f http://localhost:8080/api/health || exit 1

# Start both daemon and API server
CMD ["sh", "-c", "python3 daemon.py & python3 api_server.py"]
