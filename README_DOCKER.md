# 🚀 Dark Gene Pipeline - All-in-One Docker Image

Client 서버에 설치 가능한 통합 Docker 이미지입니다.

## 🎯 두 가지 사용 사례

### Use Case 1: Client 서버 설치 (주요 용도)
```
Client 서버 → FASTQ 분석 → Intermediate 파일 → Portal Upload → Review & Report
```
- Client가 자체 서버에서 WES 분석 수행
- Dashboard UI로 간편하게 분석 시작
- 완료 시 자동으로 Portal에 결과 업로드
- Portal에서 Review 및 Report 생성
- Dashboard에서 Report 조회 가능

### Use Case 2: Portal 중심 워크플로우
```
Portal에 FASTQ Upload → Daemon API로 Order Submit → 분석 수행 → Portal로 결과 전송
```
- Portal에서 FASTQ 업로드 및 Order 생성
- Daemon API가 파일 다운로드 및 분석 시작
- 완료 시 자동으로 결과 업로드
- Portal에서 통합 관리

## 📦 구성 요소

### 하나의 Docker Container에 포함:
- **Dashboard** (Flask, Port 5000)
  - 샘플 선택 및 분석 시작
  - 실시간 진행 상황 모니터링
  - Report Viewer
  
- **Daemon API** (Flask, Port 8080)
  - Portal과의 통신
  - Order Submit API (Use Case 2)
  - 파일 감시 및 자동 업로드
  - 상태 모니터링 API
  
- **Nextflow**
  - Multiple sample 병렬 처리 (자동)
  - Channel 기반 워크플로우

- **Supervisor**
  - 모든 프로세스 관리
  - 자동 재시작

## 🚀 빠른 설치

### 1. 자동 설치 (권장)
```bash
git clone <repository>
cd dark_gene_pipeline
./install.sh
```

설치 스크립트가 대화형으로 다음을 설정합니다:
- Data 디렉토리
- Reference 디렉토리
- Portal API 정보
- Docker 이미지 빌드 및 실행

### 2. 수동 설치

#### 2.1 환경 변수 설정
```bash
cp .env.example .env
nano .env
```

```.env
# Portal API Configuration
PORTAL_API_URL=https://your-portal.com/api
PORTAL_API_KEY=your_api_key_here
INSTITUTION_ID=your_institution_id

# Data directories
DATA_DIR=./data
REF_DIR=./references

# Pipeline configuration
BASE_DIR=/data
DAEMON_PORT=8080
NXF_OPTS=-Xms1g -Xmx4g
```

#### 2.2 디렉토리 생성
```bash
mkdir -p data/{fastq,analysis,output,log}
mkdir -p references
```

#### 2.3 Docker 이미지 빌드
```bash
cd docker
docker-compose build
```

#### 2.4 서비스 시작
```bash
docker-compose up -d
```

## 📊 사용법

### Use Case 1: Dashboard를 통한 분석

#### 1. FASTQ 파일 업로드
```bash
# 구조: data/fastq/<YYMM>/<sample_name>/
mkdir -p data/fastq/2601/Sample_A10
cp path/to/*_R1_*.fastq.gz data/fastq/2601/Sample_A10/
cp path/to/*_R2_*.fastq.gz data/fastq/2601/Sample_A10/
```

#### 2. Dashboard에서 분석 시작
```
브라우저에서 http://localhost:5000 접속
→ 샘플 선택
→ "Analyze" 버튼 클릭
→ "Monitor" 탭에서 진행 상황 확인
```

#### 3. 완료 후 Report 조회
```
샘플 목록에서 "📊 Report" 버튼 클릭
→ HTML Report 표시
```

#### 4. 자동 Portal 업로드
- 분석 완료 시 `analysis.completed` 마커 생성
- Daemon이 자동으로 감지하여 Portal에 결과 업로드

### Use Case 2: Portal API를 통한 Order Submit

#### 1. Order 제출
```bash
curl -X POST http://localhost:8080/api/orders/submit \
  -H "Content-Type: application/json" \
  -d '{
    "order_id": "ORD-2026-001",
    "work_dir": "2601",
    "sample_name": "Sample_B20",
    "fastq_r1_url": "https://portal.com/files/xxx_R1.fastq.gz",
    "fastq_r2_url": "https://portal.com/files/xxx_R2.fastq.gz",
    "priority": "normal",
    "notify_email": "user@example.com"
  }'
```

#### 2. 상태 확인
```bash
# Summary
curl http://localhost:8080/api/summary

# 전체 Order 목록
curl http://localhost:8080/api/orders

# 특정 상태의 Order만
curl http://localhost:8080/api/orders?status=RUNNING
```

#### 3. 로그 조회
```bash
# 전체 로그
curl http://localhost:8080/api/logs/2601/Sample_B20

# 마지막 100줄만
curl http://localhost:8080/api/logs/2601/Sample_B20?tail=100

# 실시간 스트리밍 (SSE)
curl http://localhost:8080/api/logs/2601/Sample_B20?follow=true
```

## 🔌 API 엔드포인트

### Dashboard API (Port 5000)

| Method | Endpoint | 설명 |
|--------|----------|------|
| GET | `/` | Dashboard UI |
| GET | `/health` | Health check |
| POST | `/start` | 분석 시작 |
| POST | `/stop` | 분석 중지 |
| GET | `/status` | 분석 상태 |
| GET | `/view_report/<work_dir>/<sample_name>` | Report 조회 |
| GET | `/api/report/<work_dir>/<sample_name>` | Report 데이터 (JSON) |
| GET | `/download/<filepath>` | 파일 다운로드 |

### Daemon API (Port 8080)

| Method | Endpoint | 설명 |
|--------|----------|------|
| GET | `/api/health` | Health check |
| GET | `/api/summary` | 분석 요약 (Portal용) |
| GET | `/api/status` | 전체 상태 |
| GET | `/api/orders` | Order 목록 |
| POST | `/api/orders/submit` | Order 제출 (Use Case 2) |
| POST | `/api/orders/<id>/retry` | Order 재시도 |
| GET | `/api/logs/<work_dir>/<sample_name>` | 로그 조회 |
| GET | `/api/logs/<work_dir>/<sample_name>?follow=true` | 실시간 로그 스트리밍 |
| GET | `/api/logs/list` | 모든 로그 목록 |

## 📁 디렉토리 구조

```
/data/                              # 마운트된 데이터 볼륨
├── fastq/                          # 입력 FASTQ
│   └── <YYMM>/                    # Work directory (예: 2601)
│       └── <sample_name>/         # 샘플별 디렉토리
│           ├── *_R1_*.fastq.gz   # Read 1
│           ├── *_R2_*.fastq.gz   # Read 2
│           └── analysis.completed # 완료 마커
│
├── analysis/                       # 중간 파일 (Intermediate)
│   └── <YYMM>/
│       └── <sample_name>/
│           ├── alignment/
│           ├── variant/
│           ├── cnv/
│           └── ...
│
├── output/                         # 최종 결과 (Portal 전송용)
│   └── <YYMM>/
│       └── <sample_name>/
│           ├── final_report.html
│           ├── variants.vcf.gz
│           └── ...
│
└── log/                           # Nextflow 로그
    └── <YYMM>/
        └── <sample_name>/
            └── .nextflow.log
```

## 🛠️ 관리 명령어

### 서비스 관리
```bash
# 상태 확인
docker ps | grep dark-gene-pipeline

# 로그 확인
docker-compose -f docker/docker-compose.yml logs -f

# Dashboard 로그만
docker-compose -f docker/docker-compose.yml logs -f dashboard

# Daemon 로그만
docker-compose -f docker/docker-compose.yml logs -f daemon

# 재시작
docker-compose -f docker/docker-compose.yml restart

# 중지
docker-compose -f docker/docker-compose.yml stop

# 완전 제거
docker-compose -f docker/docker-compose.yml down
```

### 디버깅
```bash
# 컨테이너 내부 접속
docker exec -it dark-gene-pipeline /bin/bash

# Supervisor 상태 확인
docker exec dark-gene-pipeline supervisorctl status

# 특정 프로세스 재시작
docker exec dark-gene-pipeline supervisorctl restart dashboard
docker exec dark-gene-pipeline supervisorctl restart daemon
```

### 업데이트
```bash
# 코드 업데이트 후
docker-compose -f docker/docker-compose.yml down
docker-compose -f docker/docker-compose.yml build --no-cache
docker-compose -f docker/docker-compose.yml up -d
```

## 🔐 보안 고려사항

1. **API Key 관리**
   - `.env` 파일은 `.gitignore`에 추가
   - Production 환경에서는 Secret Manager 사용 권장

2. **네트워크 설정**
   - Dashboard(5000)와 Daemon(8080) 포트를 방화벽에서 적절히 제한
   - HTTPS reverse proxy (nginx) 사용 권장

3. **데이터 보안**
   - FASTQ 파일 및 결과 파일에 적절한 권한 설정
   - 정기적인 백업

## 📊 모니터링

### Health Check
```bash
# 전체 Health Check
curl http://localhost:8080/api/health
curl http://localhost:5000/health

# Summary 확인
curl http://localhost:8080/api/summary | jq .
```

### 리소스 사용량
```bash
# Docker 컨테이너 리소스
docker stats dark-gene-pipeline

# 디스크 사용량
du -sh data/*
```

## 🐛 문제 해결

### 컨테이너가 시작되지 않음
```bash
# 로그 확인
docker logs dark-gene-pipeline

# Supervisor 로그 확인
docker exec dark-gene-pipeline cat /var/log/supervisor/supervisord.log
```

### Nextflow 에러
```bash
# Nextflow 로그 확인
cat data/log/<work_dir>/<sample_name>/.nextflow.log

# API로 확인
curl http://localhost:8080/api/logs/<work_dir>/<sample_name>
```

### Portal 연결 실패
```bash
# 환경 변수 확인
docker exec dark-gene-pipeline env | grep PORTAL

# 네트워크 테스트
docker exec dark-gene-pipeline curl -v https://your-portal.com/api
```

## 📚 추가 문서

- `BUILD_GUIDE.md` - 상세 빌드 가이드
- `API_REFERENCE.md` - API 전체 레퍼런스
- `TROUBLESHOOTING.md` - 문제 해결 가이드
- `SYSTEM_STATUS.md` - 시스템 상태 보고서

## 🤝 지원

문제가 발생하면:
1. 로그 확인: `docker-compose logs -f`
2. Health check: `curl http://localhost:8080/api/health`
3. 문서 참조: 위의 추가 문서 섹션

---

**버전**: 1.0.0  
**최종 업데이트**: 2026-01-21
