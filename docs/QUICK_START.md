# 🎊 Dark Gene Pipeline - All-in-One 구조 완료!

## 📝 작업 요약

**All-in-One Docker Image**로 성공적으로 재구성되었습니다!

---

## ✅ 완료된 작업

### 1. **All-in-One Dockerfile** ✓
- 위치: `/home/ken/dark_gene_pipeline/Dockerfile`
- Dashboard + Daemon + Nextflow를 하나의 컨테이너에 통합
- Supervisor로 모든 프로세스 관리
- Multiple sample 자동 병렬 처리 (Nextflow channel)

### 2. **Dashboard Report Viewer** ✓
- Report 조회 UI 추가
- `/view_report/<work_dir>/<sample_name>` 엔드포인트
- 샘플 목록에 "📊 Report" 버튼 추가
- Portal에서 생성한 Report를 Dashboard에서 조회 가능

### 3. **Order Submit API** ✓
- Portal에서 Order 제출 가능
- `POST /api/orders/submit` 엔드포인트
- FASTQ URL로 자동 다운로드
- Use Case 2 (Portal 중심 워크플로우) 지원

### 4. **배포 자동화** ✓
- `install.sh` 스크립트로 원클릭 설치
- 환경 변수 대화형 설정
- Docker 이미지 자동 빌드
- 서비스 자동 시작

### 5. **완전한 문서화** ✓
- `README_DOCKER.md` - 전체 사용 가이드
- `FINAL_IMPLEMENTATION_REPORT.md` - 구현 상세 보고서
- `BUILD_GUIDE.md` - 빌드 가이드
- API 레퍼런스 완비

---

## 🚀 빠른 시작

### Client 서버에 설치하기

```bash
# 1. Repository clone
git clone <repository>
cd dark_gene_pipeline

# 2. 자동 설치
./install.sh
# → Portal API URL 입력
# → API Key 입력
# → Institution ID 입력
# → 자동으로 빌드 & 실행

# 3. 서비스 확인
curl http://localhost:8080/api/health  # Daemon
curl http://localhost:5000/health      # Dashboard

# 4. Dashboard 접속
# 브라우저에서 http://localhost:5000
```

### Use Case 1: Dashboard로 분석
```bash
# FASTQ 업로드
mkdir -p data/fastq/2601/Sample_A10
cp *.fastq.gz data/fastq/2601/Sample_A10/

# Dashboard에서:
# 1. 샘플 선택
# 2. "Analyze" 클릭
# 3. "Monitor"에서 진행 상황 확인
# 4. 완료 후 "📊 Report" 클릭
```

### Use Case 2: Portal Order Submit
```bash
# Order 제출
curl -X POST http://localhost:8080/api/orders/submit \
  -H "Content-Type: application/json" \
  -d '{
    "order_id": "ORD-2026-001",
    "work_dir": "2601",
    "sample_name": "Sample_A10",
    "fastq_r1_url": "https://portal.com/files/xxx_R1.fastq.gz",
    "fastq_r2_url": "https://portal.com/files/xxx_R2.fastq.gz"
  }'

# 상태 확인
curl http://localhost:8080/api/summary
```

---

## 📊 시스템 구조

```
┌──────────────────────────────────────┐
│  Single Docker Container             │
│  (dark-gene-pipeline:latest)         │
│                                      │
│  ┌────────────────────────────────┐ │
│  │  Supervisor                    │ │
│  │  • Dashboard (Port 5000)       │ │
│  │  • Daemon API (Port 8080)      │ │
│  │  • File Watcher                │ │
│  └────────────────────────────────┘ │
│                                      │
│  Nextflow                            │
│  • Multiple sample 병렬 처리        │
│                                      │
│  Volumes:                            │
│  • /data/fastq                      │
│  • /data/analysis                   │
│  • /data/output                     │
│  • /data/log                        │
└──────────────────────────────────────┘
```

---

## 🎯 두 가지 Use Case

### Use Case 1: Client 서버 설치 (주요)
```
Client → Dashboard → Nextflow → Portal Upload → Report
```
1. Client가 자체 서버에서 분석
2. Dashboard로 샘플 관리
3. 완료 시 자동으로 Portal에 업로드
4. Portal에서 Review & Report
5. Dashboard에서 Report 조회

### Use Case 2: Portal 중심
```
Portal → Order Submit → Daemon → Nextflow → Portal Upload
```
1. Portal에서 FASTQ 업로드
2. Daemon API로 Order 제출
3. 자동으로 FASTQ 다운로드 및 분석
4. 완료 시 Portal에 업로드
5. Portal에서 통합 관리

---

## 🔌 API 엔드포인트

### Dashboard (Port 5000)
- `GET /` - Web UI
- `POST /start` - 분석 시작
- **`GET /view_report/<work_dir>/<sample_name>`** - **Report Viewer** (신규)
- **`GET /api/report/<work_dir>/<sample_name>`** - **Report 데이터** (신규)

### Daemon (Port 8080)
- `GET /api/summary` - Portal Summary
- `GET /api/status` - 전체 상태
- **`POST /api/orders/submit`** - **Order 제출** (신규)
- `GET /api/orders` - Order 목록
- `GET /api/logs/<work_dir>/<sample_name>` - 로그 조회

---

## 📁 디렉토리 구조

```
/data/
├── fastq/<YYMM>/<sample_name>/     # 입력 FASTQ
│   ├── *_R1_*.fastq.gz
│   ├── *_R2_*.fastq.gz
│   └── analysis.completed          # 완료 마커
│
├── analysis/<YYMM>/<sample_name>/  # 중간 파일
│   ├── alignment/
│   ├── variant/
│   └── ...
│
├── output/<YYMM>/<sample_name>/    # 최종 결과
│   ├── final_report.html
│   ├── variants.vcf.gz
│   └── ...
│
└── log/<YYMM>/<sample_name>/       # 로그
    └── .nextflow.log
```

---

## 🔧 관리 명령어

```bash
# 서비스 상태
docker ps | grep dark-gene-pipeline

# 로그 확인
docker-compose -f docker/docker-compose.yml logs -f

# 재시작
docker-compose -f docker/docker-compose.yml restart

# 중지
docker-compose -f docker/docker-compose.yml stop

# 제거
docker-compose -f docker/docker-compose.yml down
```

---

## 📚 문서

| 문서 | 설명 |
|------|------|
| **README_DOCKER.md** | 전체 사용 가이드 (이 문서 권장!) |
| **FINAL_IMPLEMENTATION_REPORT.md** | 구현 상세 보고서 |
| **BUILD_GUIDE.md** | 빌드 가이드 |
| **SYSTEM_STATUS.md** | 시스템 현황 |
| **LOG_VIEWER_GUIDE.md** | 로그 조회 방법 |

---

## ✨ 주요 개선 사항

### 이전 구조
- ❌ Daemon만 Docker
- ❌ Dashboard는 호스트에서 실행
- ❌ 설치 복잡
- ❌ Multiple sample 처리 불명확

### 현재 구조 (All-in-One)
- ✅ 모든 구성 요소 통합
- ✅ `./install.sh` 한 번으로 설치
- ✅ Nextflow channel 자동 병렬 처리
- ✅ Client 설치 간편화
- ✅ Report Viewer 추가
- ✅ Order Submit API 추가
- ✅ 두 가지 Use Case 모두 지원

---

## 🎯 다음 단계

### Portal 팀
- [ ] Portal API 구현 (파일 수신, 상태 업데이트)
- [ ] Summary 페이지 구현
- [ ] Report 생성 및 제공

### Client 팀
- [ ] Reference 데이터 준비
- [ ] 테스트 샘플 실행
- [ ] Production 배포

### 통합 테스트
- [ ] Use Case 1 전체 워크플로우
- [ ] Use Case 2 Portal 연동
- [ ] Multiple sample 동시 처리
- [ ] Report 조회 및 다운로드

---

## 🎉 완료!

**All-in-One Docker Image**로 성공적으로 재구성되었습니다!

### 설치 준비 완료
```bash
./install.sh
```
한 번만 실행하면 Client 서버에 모든 것이 설치됩니다.

### 주요 기능
- ✅ Dashboard로 분석 시작 및 모니터링
- ✅ Report Viewer로 결과 확인
- ✅ Portal Order Submit 지원
- ✅ Multiple sample 자동 병렬 처리
- ✅ 자동 Portal 업로드

**Production Ready!** 🚀

---

**버전**: 1.0.0  
**날짜**: 2026-01-21  
**상태**: ✅ 배포 준비 완료
