# Dark Gene Pipeline - 전체 시스템 구조

## 🏗️ 아키텍처 개요

```
┌─────────────────────────────────────────────────────────────┐
│                         Portal                              │
│                  (portal.genolyx.com)                       │
│  - Analysis Summary Dashboard (첨부 이미지와 같은 UI)        │
│  - 결과 파일 조회 및 다운로드                                │
│  - 사용자 관리                                              │
└─────────────────┬───────────────────────────────────────────┘
                  │ REST API
                  │ - 상태 조회
                  │ - 파일 업로드
                  │ - 진행률 보고
                  ▼
┌─────────────────────────────────────────────────────────────┐
│              Portal Daemon (Docker Container)               │
│  - 분석 완료 감지 (File System Watcher)                      │
│  - 자동 결과 업로드                                          │
│  - 실시간 상태 모니터링                                      │
│  - API Server (localhost:8080)                             │
└─────────────────┬───────────────────────────────────────────┘
                  │ File System Access (Read-only)
                  │
┌─────────────────▼───────────────────────────────────────────┐
│            Dark Gene Pipeline (Host System)                 │
│                                                             │
│  ┌─────────────────────────────────────────────────────┐   │
│  │            Dashboard (Flask App)                     │   │
│  │  - 샘플 선택 및 분석 시작                             │   │
│  │  - 실시간 진행 모니터링                               │   │
│  │  - 로컬 결과 확인                                    │   │
│  │  localhost:5000                                     │   │
│  └─────────────────┬───────────────────────────────────┘   │
│                    │                                         │
│  ┌─────────────────▼───────────────────────────────────┐   │
│  │          Nextflow Pipeline                          │   │
│  │  - BWA, GATK, Manta, Paraphase                      │   │
│  │  - Docker/Singularity containers                    │   │
│  │  - Micromamba 환경 관리                              │   │
│  └─────────────────┬───────────────────────────────────┘   │
│                    │                                         │
│  ┌─────────────────▼───────────────────────────────────┐   │
│  │         Directory Structure                         │   │
│  │  fastq/<work_dir>/<sample>/     (원본)              │   │
│  │  analysis/<work_dir>/<sample>/  (intermediate)      │   │
│  │  output/<work_dir>/<sample>/    (final results)     │   │
│  │  log/<work_dir>/<sample>/       (logs)              │   │
│  │  work/                          (Nextflow temp)     │   │
│  └─────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

---

## 🐳 Docker 구성

### 1. Dark Gene Pipeline (현재 구조)
- **타입**: Host 시스템에서 직접 실행
- **Docker 사용**: 각 프로세스별 container 사용
  - `quay.io/biocontainers/bwa:0.7.17`
  - `broadinstitute/gatk:4.4.0.0`
  - `quay.io/biocontainers/manta:1.6.0`
  - `quay.io/biocontainers/paraphase:3.1.0`
  - etc.
- **실행 방식**: Nextflow가 필요에 따라 Docker pull & run
- **프로파일**: `docker` 또는 `singularity`

### 2. Portal Daemon (새로 추가)
- **타입**: 독립적인 Docker container
- **Base Image**: `python:3.11-slim`
- **실행**: Docker Compose 또는 단독 실행
- **Volume Mount**: 파이프라인 디렉토리 (read-only)
- **네트워크**: Host 네트워크 또는 브리지

---

## 🔄 데이터 흐름

### 1. 분석 시작 (Dashboard → Pipeline)
```
사용자 (Dashboard)
   ↓ 샘플 선택 & Start
Dashboard (Flask)
   ↓ Nextflow 실행
Pipeline (Nextflow + Docker containers)
   ↓ 분석 수행
Results → output/<work_dir>/<sample>/
   ↓ analysis.completed 생성
fastq/<work_dir>/<sample>/analysis.completed
```

### 2. 결과 업로드 (Daemon → Portal)
```
Daemon (File Watcher)
   ↓ analysis.completed 감지
Daemon (Scanner)
   ↓ output/ 디렉토리 스캔
Daemon (Uploader)
   ↓ Portal API 호출
Portal
   ↓ 파일 저장 & DB 기록
Portal Dashboard
   ↓ 사용자에게 표시
```

### 3. 상태 모니터링 (Portal ← Daemon)
```
Portal (Frontend)
   ↓ GET /api/analysis/summary
Portal (Backend)
   ↓ 내부 또는 Daemon API 호출
Daemon API (GET /api/summary)
   ↓ 디렉토리 스캔 & 상태 집계
Response: {
    today_requested: 0,
    running: 0,
    today_completed: 0,
    ...
}
```

---

## 📁 파일 시스템 구조

### Pipeline 관점
```
/home/ken/dark_gene_pipeline/
├── fastq/2601/Sample_A10/
│   ├── *_R1_*.fq.gz
│   ├── *_R2_*.fq.gz
│   └── analysis.completed ← Daemon 감시
├── analysis/2601/Sample_A10/
│   ├── alignment/*.bam
│   ├── variant/*.vcf.gz
│   └── pipeline_info/trace.txt ← Daemon 파싱 (진행률)
├── output/2601/Sample_A10/ ← Daemon 업로드 대상
│   ├── summary/*.txt
│   ├── snapshots/*.html
│   ├── vcf/*.vcf.gz
│   └── bam/*.bam
└── log/2601/Sample_A10/
    └── nextflow.log
```

### Daemon 관점 (Docker Volume Mount)
```
/pipeline/ (inside container)
├── fastq/     → /home/ken/.../fastq:ro
├── analysis/  → /home/ken/.../analysis:ro
├── output/    → /home/ken/.../output:ro
└── log/       → /home/ken/.../log:ro
```

---

## 🚀 배포 시나리오

### 시나리오 1: All-in-One (현재 권장)
```
Server A (분석 서버)
├── Dark Gene Pipeline (Host)
│   ├── Nextflow
│   ├── Docker containers (per process)
│   └── Dashboard (Flask, port 5000)
└── Portal Daemon (Docker, port 8080)
    └── Volume mounts → Pipeline directories
```

**장점:**
- 단일 서버 관리
- 파일 시스템 직접 접근 (빠름)
- 네트워크 지연 없음

**단점:**
- 서버 리소스 경합 가능
- 확장성 제한

---

### 시나리오 2: 분리 배포 (대규모 환경)
```
Server A (분석 서버)
└── Dark Gene Pipeline (Host)
    ├── Nextflow
    └── Dashboard (port 5000)
    
Server B (Daemon 서버)
└── Portal Daemon (Docker, port 8080)
    └── NFS/SMB mount → Server A directories

Server C (Portal 서버)
└── Portal Application
    └── API calls → Daemon (Server B)
```

**장점:**
- 역할 분리 (분석 vs 업로드)
- 확장 가능
- Portal 독립적 운영

**단점:**
- 네트워크 파일 시스템 필요
- 관리 복잡도 증가

---

## 🔌 Portal API 스펙

### Portal이 제공해야 하는 API

#### 1. 주문 생성
```http
POST /api/v1/analysis/order
Authorization: Bearer {api_key}

Request:
{
    "work_dir": "2601",
    "sample_name": "Sample_A10",
    "status": "WAITING",
    "institution_id": "genolyx"
}

Response:
{
    "order_id": "ORD-2601-001",
    "status": "created"
}
```

#### 2. 상태 업데이트
```http
PUT /api/v1/analysis/order/{order_id}
Authorization: Bearer {api_key}

Request:
{
    "status": "RUNNING",
    "progress": 45,
    "details": "CALL_VARIANTS in progress"
}

Response:
{
    "status": "updated"
}
```

#### 3. 파일 업로드
```http
POST /api/v1/analysis/upload
Authorization: Bearer {api_key}
Content-Type: multipart/form-data

Form Data:
- file: (binary)
- order_id: ORD-2601-001
- file_type: summary
- file_size: 12345

Response:
{
    "file_id": "FILE-001",
    "url": "https://portal.../files/FILE-001"
}
```

#### 4. 완료 보고
```http
POST /api/v1/analysis/order/{order_id}/complete
Authorization: Bearer {api_key}

Request:
{
    "status": "COMPLETED",
    "duration": "2h 15m",
    "completed_at": "2026-01-21T14:30:00Z"
}

Response:
{
    "status": "completed"
}
```

---

## 📊 상태 동기화

### Portal Summary 페이지 (첨부 이미지)

**데이터 소스 옵션:**

#### 옵션 1: Daemon API 직접 호출
```javascript
// Portal Frontend
async function updateSummary() {
    const response = await fetch('http://daemon-server:8080/api/summary');
    const data = await response.json();
    
    // UI 업데이트
    document.getElementById('today-requested').textContent = data.today_requested;
    document.getElementById('running').textContent = `${data.running}/8`;
    // ...
}
```

#### 옵션 2: Portal DB 기반 (권장)
```
Daemon → Portal API → Portal DB
Portal Frontend → Portal Backend → Portal DB
```

**장점:**
- Portal이 데이터 소유
- 히스토리 관리 용이
- 복수 Daemon 지원 가능

---

## 🔧 설정 예시

### 1. Daemon 환경 변수 (`.env`)
```bash
PORTAL_URL=https://portal.genolyx.com
PORTAL_API_KEY=pk_live_abc123xyz789
INSTITUTION_ID=genolyx
PIPELINE_BASE_DIR=/pipeline
```

### 2. Docker Compose 실행
```bash
cd /home/ken/dark_gene_pipeline/daemon

# .env 파일 생성
cp .env.example .env
vim .env  # API Key 설정

# 실행
docker-compose up -d

# 확인
docker-compose ps
curl http://localhost:8080/api/health
```

### 3. Dashboard 실행 (별도)
```bash
cd /home/ken/dark_gene_pipeline/dashboard
python3 app.py
# http://localhost:5000
```

---

## 🎯 통합 테스트 시나리오

### 1. 전체 흐름 테스트
```bash
# 1. Daemon 시작
cd daemon && docker-compose up -d

# 2. Dashboard 시작
cd dashboard && python3 app.py &

# 3. 샘플 분석 시작 (Dashboard)
# http://localhost:5000 → Select Sample → Analyze

# 4. Daemon 로그 모니터링
docker-compose logs -f daemon

# 5. Daemon API 확인
curl http://localhost:8080/api/summary

# 6. Portal 확인 (실제 Portal 필요)
# https://portal.genolyx.com/analysis/summary
```

### 2. 수동 업로드 테스트
```bash
# 분석 완료 마커 생성 (수동)
touch /home/ken/dark_gene_pipeline/fastq/2601/Sample_Test/analysis.completed

# Daemon이 자동 감지 및 업로드 (약 5초 내)
docker-compose logs -f daemon | grep "Sample_Test"
```

---

## 📝 요약

### Dark Gene Pipeline
- ✅ Host 시스템에서 실행
- ✅ Nextflow + Docker containers (per process)
- ✅ Dashboard (Flask, localhost:5000)

### Portal Daemon
- ✅ **별도 Docker container**
- ✅ 독립적 이미지 (`dark-gene-daemon:latest`)
- ✅ Pipeline 디렉토리 마운트 (read-only)
- ✅ Portal API 클라이언트
- ✅ REST API 서버 (localhost:8080)

### 통합 포인트
- ✅ File System: `analysis.completed` 마커
- ✅ API: Portal REST API
- ✅ 네트워크: Docker bridge 또는 host

**이제 Portal API만 준비되면 전체 시스템이 작동합니다!** 🚀
