# Dark Gene Pipeline - Portal Daemon

Portal과 연동하여 분석 결과를 자동으로 업로드하고 상태를 모니터링하는 Daemon 서비스입니다.

## 🎯 주요 기능

### 1. 자동 결과 업로드
- 분석 완료된 샘플 자동 감지 (`analysis.completed` 마커)
- Portal API를 통해 결과 파일 자동 업로드
- Summary, Snapshots, VCF, BAM 등 모든 결과 포함

### 2. 실시간 상태 모니터링
- 실행 중인 샘플 진행률 추적
- Portal에 주기적으로 상태 보고
- Nextflow trace 파일 기반 진행률 계산

### 3. Portal API 제공
```
GET  /api/summary          # 분석 상태 요약 (Portal Summary 페이지)
GET  /api/status           # 전체 상태 조회
GET  /api/orders           # 주문 목록
POST /api/orders/{id}/retry # 재시도
POST /api/samples/scan     # 수동 스캔
```

### 4. File System Watching
- `fastq/` 디렉토리 실시간 감시
- `analysis.completed` 생성 즉시 업로드 시작

---

## 🚀 설치 및 실행

### 방법 1: Docker Compose (권장)

#### 1. 환경 변수 설정
```bash
cd daemon
cp .env.example .env

# .env 파일 편집
vim .env
```

```env
PORTAL_URL=https://portal.genolyx.com
PORTAL_API_KEY=your_actual_api_key_here
INSTITUTION_ID=genolyx
```

#### 2. Docker 이미지 빌드 및 실행
```bash
# 빌드
docker-compose build

# 실행
docker-compose up -d

# 로그 확인
docker-compose logs -f daemon
```

#### 3. 상태 확인
```bash
# 헬스 체크
curl http://localhost:8080/api/health

# 분석 요약
curl http://localhost:8080/api/summary
```

---

### 방법 2: Docker 단독 실행

```bash
# 이미지 빌드
docker build -t dark-gene-daemon:latest .

# 실행
docker run -d \
  --name dark-gene-daemon \
  --restart unless-stopped \
  -v /home/ken/dark_gene_pipeline/fastq:/pipeline/fastq:ro \
  -v /home/ken/dark_gene_pipeline/analysis:/pipeline/analysis:ro \
  -v /home/ken/dark_gene_pipeline/output:/pipeline/output:ro \
  -v /home/ken/dark_gene_pipeline/log:/pipeline/log:ro \
  -e PORTAL_URL=https://portal.genolyx.com \
  -e PORTAL_API_KEY=your_api_key_here \
  -e INSTITUTION_ID=genolyx \
  -p 8080:8080 \
  dark-gene-daemon:latest

# 로그 확인
docker logs -f dark-gene-daemon
```

---

### 방법 3: 로컬 실행 (개발용)

```bash
# 의존성 설치
pip install -r requirements.txt

# 환경 변수 설정
export PIPELINE_BASE_DIR=/home/ken/dark_gene_pipeline
export PORTAL_URL=https://portal.genolyx.com
export PORTAL_API_KEY=your_api_key_here
export INSTITUTION_ID=genolyx

# Daemon 실행 (백그라운드)
python3 daemon.py &

# API Server 실행
python3 api_server.py
```

---

## 📡 Portal API 연동

### Portal에서 요구하는 API 엔드포인트

Daemon이 호출하는 Portal API:

```python
# 1. 분석 주문 생성
POST /api/v1/analysis/order
{
    "work_dir": "2601",
    "sample_name": "Sample_A10",
    "status": "WAITING",
    "institution_id": "genolyx"
}

# 2. 상태 업데이트
PUT /api/v1/analysis/order/{order_id}
{
    "status": "RUNNING",
    "progress": 45,
    "updated_at": "2026-01-21T..."
}

# 3. 파일 업로드
POST /api/v1/analysis/upload
- multipart/form-data
- file: 실제 파일
- order_id: 주문 ID
- file_type: summary, vcf, bam, etc.

# 4. 완료 보고
POST /api/v1/analysis/order/{order_id}/complete
{
    "status": "COMPLETED",
    "duration": "2h 15m",
    "completed_at": "2026-01-21T..."
}
```

---

## 📊 분석 상태 Summary API

Portal Summary 페이지에서 호출:

```bash
GET http://localhost:8080/api/summary
```

**Response:**
```json
{
    "today_requested": 0,
    "running": 0,
    "queue_waiting": 0,
    "today_completed": 0,
    "today_failed": 0,
    "total_requested": 72,
    "total_completed": 72,
    "total_failed": 0
}
```

Portal에서 이 데이터를 받아 첨부 이미지와 같은 UI 구성:
- Today Requested: `today_requested`
- Running: `running` / 8
- Queue Waiting: `queue_waiting`
- Today Completed: `today_completed`
- Today Failed: `today_failed`
- Total Requested: `total_requested`
- Total Completed: `total_completed`
- Total Failed: `total_failed`

---

## 🔄 동작 원리

### 1. 샘플 스캔
```
매 30초마다:
1. fastq/ 디렉토리 스캔
2. analysis.completed 마커 확인
3. 완료된 샘플 → 업로드 큐에 추가
```

### 2. 진행 중 샘플 모니터링
```
매 10초마다:
1. analysis/ 디렉토리 확인
2. trace.txt 파일 파싱
3. 진행률 계산 (0-100%)
4. Portal에 상태 보고
```

### 3. 실시간 감지
```
File System Watcher:
1. fastq/를 실시간 감시
2. analysis.completed 생성 감지
3. 즉시 업로드 프로세스 시작
```

### 4. 업로드 프로세스
```
1. Portal에 order_id 생성 (없으면)
2. output/ 디렉토리 스캔
3. 파일 타입별로 업로드
   - summary/
   - snapshots/
   - vcf/
   - bam/
   - cnv/
   - sv/
   - repeat/
   - pseudogene/
4. 완료 보고
5. state.json에 상태 저장
```

---

## 🛠️ 관리 명령어

### Docker Compose

```bash
# 시작
docker-compose up -d

# 중지
docker-compose stop

# 재시작
docker-compose restart

# 로그 확인
docker-compose logs -f

# 상태 확인
docker-compose ps

# 완전 삭제 (볼륨 포함)
docker-compose down -v
```

### 수동 제어

```bash
# 수동 스캔 트리거
curl -X POST http://localhost:8080/api/samples/scan

# 특정 주문 재시도
curl -X POST http://localhost:8080/api/orders/{order_id}/retry

# 상태 조회
curl http://localhost:8080/api/status
```

---

## 📝 상태 파일 구조

`/var/lib/daemon/state.json`:
```json
{
    "orders": {
        "2601/Sample_A10": {
            "order_id": "ORD-2601-001",
            "created_at": "2026-01-21T10:00:00",
            "uploaded": true,
            "uploaded_at": "2026-01-21T12:15:00"
        }
    }
}
```

---

## 🔍 트러블슈팅

### 1. Daemon이 Portal에 연결 못함
```bash
# 로그 확인
docker logs dark-gene-daemon | grep "ERROR"

# API Key 확인
docker exec dark-gene-daemon env | grep PORTAL_API_KEY

# Portal 연결 테스트
docker exec dark-gene-daemon curl -H "Authorization: Bearer $API_KEY" $PORTAL_URL/api/v1/health
```

### 2. 업로드 실패
```bash
# 상태 파일 확인
docker exec dark-gene-daemon cat /var/lib/daemon/state.json

# 수동 스캔
curl -X POST http://localhost:8080/api/samples/scan
```

### 3. 진행률이 업데이트 안됨
```bash
# trace 파일 확인
ls -l /home/ken/dark_gene_pipeline/analysis/*/*/pipeline_info/trace.txt

# Daemon 재시작
docker-compose restart
```

---

## 🔐 보안

### API Key 관리
- `.env` 파일에 저장 (git에 커밋 금지)
- Docker secrets 사용 권장 (프로덕션)
- 주기적 로테이션

### 네트워크
- 내부 네트워크만 접근 가능하도록 설정
- Portal만 Daemon API 호출 가능

### 파일 권한
- 읽기 전용 마운트 (`ro` 플래그)
- Daemon은 파이프라인 파일 수정 불가

---

## 📈 모니터링

### Prometheus Metrics (향후 추가 가능)
```
dark_gene_samples_total
dark_gene_samples_running
dark_gene_samples_completed
dark_gene_upload_errors_total
```

### Grafana Dashboard
- 실시간 분석 현황
- 업로드 성공률
- 평균 분석 시간

---

## 🚀 다음 단계

1. **Portal API Key 발급**
   - Portal 관리자에게 요청
   - Institution ID 확인

2. **Daemon 배포**
   ```bash
   cd daemon
   vim .env  # API Key 설정
   docker-compose up -d
   ```

3. **Portal에서 확인**
   - Summary 페이지 접속
   - 분석 상태 확인

4. **테스트**
   - 샘플 분석 실행
   - 완료 후 Portal에서 결과 확인
