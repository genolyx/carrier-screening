#!/bin/bash
#
# Dark Gene Pipeline - Quick Deploy Script
# Client 서버에 설치하기 위한 스크립트
#

set -e

echo "======================================"
echo "Dark Gene Pipeline - All-in-One Setup"
echo "======================================"
echo ""

# Check Docker
if ! command -v docker &> /dev/null; then
    echo "❌ Docker is not installed. Please install Docker first."
    exit 1
fi

if ! command -v docker-compose &> /dev/null; then
    echo "❌ Docker Compose is not installed. Please install Docker Compose first."
    exit 1
fi

echo "✅ Docker and Docker Compose are installed"
echo ""

# Configuration
echo "📋 Configuration"
echo "----------------"

read -p "Data directory [./data]: " DATA_DIR
DATA_DIR=${DATA_DIR:-./data}

read -p "Reference directory [./references]: " REF_DIR
REF_DIR=${REF_DIR:-./references}

read -p "Portal API URL [https://portal.example.com/api]: " PORTAL_API_URL
PORTAL_API_URL=${PORTAL_API_URL:-https://portal.example.com/api}

read -p "Portal API Key: " PORTAL_API_KEY

read -p "Institution ID: " INSTITUTION_ID

echo ""

# Create .env file
echo "📝 Creating .env file..."
cat > .env << EOF
# Portal API Configuration
PORTAL_API_URL=${PORTAL_API_URL}
PORTAL_API_KEY=${PORTAL_API_KEY}
INSTITUTION_ID=${INSTITUTION_ID}

# Data directories
DATA_DIR=${DATA_DIR}
REF_DIR=${REF_DIR}

# Pipeline configuration
BASE_DIR=/data
DAEMON_PORT=8080
NXF_OPTS=-Xms1g -Xmx4g

# Log level
LOG_LEVEL=INFO
EOF

echo "✅ .env file created"
echo ""

# Create data directories
echo "📁 Creating data directories..."
mkdir -p ${DATA_DIR}/{fastq,analysis,output,log}
echo "✅ Data directories created"
echo ""

# Build Docker image
echo "🐳 Building Docker image..."
docker-compose -f docker/docker-compose.yml build
echo "✅ Docker image built"
echo ""

# Start services
echo "🚀 Starting services..."
docker-compose -f docker/docker-compose.yml up -d
echo "✅ Services started"
echo ""

# Wait for services to be ready
echo "⏳ Waiting for services to be ready..."
sleep 10

# Health check
echo "🏥 Health check..."
if curl -f http://localhost:8080/api/health > /dev/null 2>&1; then
    echo "✅ Daemon API is healthy"
else
    echo "⚠️  Daemon API is not responding yet"
fi

if curl -f http://localhost:5000/health > /dev/null 2>&1; then
    echo "✅ Dashboard is healthy"
else
    echo "⚠️  Dashboard is not responding yet"
fi

echo ""
echo "======================================"
echo "🎉 Installation Complete!"
echo "======================================"
echo ""
echo "Services:"
echo "  Dashboard: http://localhost:5000"
echo "  Daemon API: http://localhost:8080"
echo ""
echo "Data directories:"
echo "  FASTQ: ${DATA_DIR}/fastq"
echo "  Analysis: ${DATA_DIR}/analysis"
echo "  Output: ${DATA_DIR}/output"
echo "  Logs: ${DATA_DIR}/log"
echo ""
echo "Management commands:"
echo "  View logs: docker-compose -f docker/docker-compose.yml logs -f"
echo "  Stop: docker-compose -f docker/docker-compose.yml stop"
echo "  Restart: docker-compose -f docker/docker-compose.yml restart"
echo "  Remove: docker-compose -f docker/docker-compose.yml down"
echo ""
echo "Next steps:"
echo "  1. Place your reference files in ${REF_DIR}/"
echo "  2. Upload FASTQ files to ${DATA_DIR}/fastq/<YYMM>/<sample_name>/"
echo "  3. Start analysis from Dashboard (http://localhost:5000)"
echo ""
