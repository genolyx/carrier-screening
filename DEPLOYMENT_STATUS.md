# 🚀 Dark Gene Pipeline Daemon - Quick Start Guide

## ✅ 빌드 및 배포 완료!

Daemon이 성공적으로 빌드되어 실행 중입니다.

## 📊 현재 상태

### 서비스 확인
```bash
# 컨테이너 상태
docker ps | grep dark-gene-daemon

# 로그 확인
docker logs -f dark-gene-daemon
```

### API 테스트
```bash
# Health Check
curl http://localhost:8080/api/health

# 전체 Status
curl http://localhost:8080/api/status | python3 -m json.tool

# Summary (Portal용)
curl http://localhost:8080/api/summary | python3 -m json.tool

# 특정 샘플의 로그 조회
curl http://localhost:8080/api/logs/2601/Sample_Test | python3 -m json.tool
```

### 웹 인터페이스
브라우저에서 접속:
- Log Viewer: http://localhost:8080/

## 🎯 다음 단계

### 1. Portal API 연동
`.env` 파일을 수정하여 실제 Portal API 정보를 입력하세요:

```bash
nano /home/ken/dark_gene_pipeline/daemon/.env
```

다음 항목을 수정:
- `PORTAL_API_URL`: Portal의 실제 API URL
- `PORTAL_API_KEY`: Portal에서 발급받은 API Key
- `INSTITUTION_ID`: 기관 ID

변경 후 서비스 재시작:
```bash
cd /home/ken/dark_gene_pipeline/daemon
docker-compose restart
```

### 2. Portal 개발
Portal에서 다음 기능을 구현하세요:

#### 2.1 Summary API 연동
```javascript
// Portal Frontend Example
fetch('http://daemon-server:8080/api/summary')
  .then(res => res.json())
  .then(data => {
    // data.today_completed, data.running 등을 화면에 표시
  });
```

#### 2.2 Status API 연동
```javascript
// 실시간 샘플 모니터링
setInterval(() => {
  fetch('http://daemon-server:8080/api/status')
    .then(res => res.json())
    .then(data => {
      // data.samples 배열을 통해 각 샘플의 상태 표시
    });
}, 5000); // 5초마다 업데이트
```

#### 2.3 Log Viewer 연동
```javascript
// 특정 샘플의 로그 조회
function viewLogs(workDir, sampleName) {
  fetch(`http://daemon-server:8080/api/logs/${workDir}/${sampleName}`)
    .then(res => res.json())
    .then(data => {
      console.log(data.log_content);
    });
}
```

#### 2.4 파일 업로드 수신 API
Daemon이 분석 완료 후 결과를 업로드할 Portal API를 구현:
```python
# Portal Backend Example (Flask)
@app.route('/api/orders/<order_id>/upload', methods=['POST'])
def upload_result(order_id):
    file = request.files['file']
    file_type = request.form['file_type']
    # 파일 저장 로직
    return jsonify({'status': 'success'})
```

### 3. 통합 테스트

#### 3.1 샘플 분석 시작
Dashboard에서 샘플 분석을 시작합니다:
```bash
# 브라우저에서 접속
http://localhost:5000
```

#### 3.2 Daemon 모니터링
Daemon이 분석 상태를 감지하는지 확인:
```bash
# Daemon 로그 확인
docker logs -f dark-gene-daemon

# Status API로 확인
watch -n 5 'curl -s http://localhost:8080/api/status | python3 -m json.tool'
```

#### 3.3 완료 후 자동 업로드 확인
- `analysis.completed` 파일 생성 시 Daemon이 자동으로 Portal에 업로드
- Portal에서 수신된 파일 확인

## 🔧 관리 명령어

### 서비스 관리
```bash
cd /home/ken/dark_gene_pipeline/daemon

# 시작
docker-compose up -d

# 중지
docker-compose stop

# 재시작
docker-compose restart

# 로그 확인
docker-compose logs -f

# 완전 종료 및 삭제
docker-compose down
```

### 코드 수정 후 재배포
```bash
cd /home/ken/dark_gene_pipeline/daemon

# 기존 컨테이너 중지 및 삭제
docker-compose down

# 이미지 재빌드 (캐시 없이)
docker-compose build --no-cache

# 재시작
docker-compose up -d

# 로그 확인
docker-compose logs -f
```

### 디버깅
```bash
# 컨테이너 내부 접속
docker exec -it dark-gene-daemon /bin/bash

# 환경 변수 확인
docker exec dark-gene-daemon env

# 파일 시스템 확인
docker exec dark-gene-daemon ls -la /pipeline/
```

## 📁 디렉토리 구조

```
/home/ken/dark_gene_pipeline/
├── daemon/                    # Daemon 소스 코드
│   ├── api_server.py         # Flask API 서버
│   ├── daemon.py             # 파일 감시 및 업로드 로직
│   ├── Dockerfile            # Docker 이미지 빌드 파일
│   ├── docker-compose.yml    # Docker Compose 설정
│   ├── requirements.txt      # Python 의존성
│   ├── .env                  # 환경 변수 (수정 필요)
│   └── BUILD_GUIDE.md        # 상세 빌드 가이드
├── fastq/                    # 입력 FASTQ 파일
│   └── <work_dir>/
│       └── <sample_name>/
│           ├── *_R1_*.fastq.gz
│           ├── *_R2_*.fastq.gz
│           └── analysis.completed  # 완료 마커
├── analysis/                 # 분석 중간 파일
│   └── <work_dir>/
│       └── <sample_name>/
├── output/                   # Portal 전송용 최종 결과
│   └── <work_dir>/
│       └── <sample_name>/
└── log/                      # Nextflow 로그
    └── <work_dir>/
        └── <sample_name>/
            └── .nextflow.log
```

## 🚨 문제 해결

### 컨테이너가 계속 재시작되는 경우
```bash
docker logs dark-gene-daemon
# 로그에서 에러 메시지 확인
```

### API 연결이 안 되는 경우
```bash
# 포트 확인
netstat -tuln | grep 8080

# 방화벽 확인
sudo ufw status
sudo ufw allow 8080
```

### Portal API 호출 실패
```bash
# .env 파일 확인
cat /home/ken/dark_gene_pipeline/daemon/.env

# 네트워크 테스트
docker exec dark-gene-daemon curl -v https://your-portal-domain.com/api
```

### 파일 감지가 안 되는 경우
```bash
# 볼륨 마운트 확인
docker inspect dark-gene-daemon | grep -A 10 Mounts

# 파일 권한 확인
ls -la /home/ken/dark_gene_pipeline/fastq/
```

## 📊 성능 모니터링

### 리소스 사용량 확인
```bash
# 컨테이너 리소스 사용량
docker stats dark-gene-daemon

# 디스크 사용량
df -h /home/ken/dark_gene_pipeline/
```

### 로그 크기 관리
```bash
# 로그 파일 크기 확인
docker inspect dark-gene-daemon | grep LogPath

# 로그 정리
docker-compose down
docker system prune -f
docker-compose up -d
```

## 🔐 보안 체크리스트

- [ ] `.env` 파일이 `.gitignore`에 포함되어 있는지 확인
- [ ] API Key가 안전하게 관리되고 있는지 확인
- [ ] Portal API가 HTTPS를 사용하는지 확인
- [ ] 컨테이너 포트(8080)가 외부에 노출되지 않도록 설정
- [ ] 정기적으로 API Key 갱신
- [ ] 로그 파일에 민감한 정보가 기록되지 않도록 확인

## 📞 지원

문제가 발생하면 다음을 확인하세요:
1. `BUILD_GUIDE.md` - 상세 빌드 가이드
2. `daemon/README.md` - Daemon 기능 설명
3. `LOG_VIEWER_GUIDE.md` - 로그 조회 방법
4. Docker 로그: `docker logs dark-gene-daemon`

---

**현재 상태**: ✅ Daemon 실행 중 (http://localhost:8080)  
**다음 단계**: Portal API Key 설정 및 통합 테스트
