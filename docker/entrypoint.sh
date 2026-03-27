#!/bin/bash
set -e

# Create necessary directories
mkdir -p /data/fastq /data/analysis /data/output /data/log

exec "$@"
