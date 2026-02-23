# Daemon Build & Deployment Guide

## 빌드 전 준비사항

### 1. 환경 설정 파일 수정
```bash
cd /home/ken/dark_gene_pipeline/daemon
nano .env
```

다음 항목들을 수정하세요:
- `PORTAL_API_URL`: Portal의 실제 API URL
- `PORTAL_API_KEY`: Portal에서 발급받은 API Key
- `INSTITUTION_ID`: 기관 ID

### 2. 디렉토리 권한 확인
```bash
# 파이프라인 디렉토리가 읽기/쓰기 가능한지 확인
ls -la /home/ken/dark_gene_pipeline/
```

## 빌드 및 배포

### Option 1: Docker Compose 사용 (권장)

#### 1단계: 이미지 빌드
```bash
cd /home/ken/dark_gene_pipeline/daemon
docker-compose build
```

#### 2단계: 서비스 시작
```bash
docker-compose up -d
```

#### 3단계: 로그 확인
```bash
docker-compose logs -f
```

#### 4단계: 서비스 중지
```bash
docker-compose down
```

### Option 2: Docker 직접 사용

#### 1단계: 이미지 빌드
```bash
cd /home/ken/dark_gene_pipeline/daemon
docker build -t dark-gene-daemon:latest .
```

#### 2단계: 컨테이너 실행
```bash
docker run -d \
  --name dark-gene-daemon \
  -p 8080:8080 \
  -v /home/ken/dark_gene_pipeline:/dark_gene_pipeline \
  --env-file .env \
  --restart unless-stopped \
  dark-gene-daemon:latest
```

#### 3단계: 로그 확인
```bash
docker logs -f dark-gene-daemon
```

#### 4단계: 컨테이너 중지/삭제
```bash
docker stop dark-gene-daemon
docker rm dark-gene-daemon
```

## 서비스 확인

### 1. Daemon API 테스트
```bash
# Health check
curl http://localhost:8080/health

# Status 확인
curl http://localhost:8080/api/status

# Summary 확인
curl http://localhost:8080/api/summary
```

### 2. Log Viewer 테스트
브라우저에서 접속:
```
http://localhost:8080/
```

### 3. 특정 샘플 로그 확인
```bash
curl http://localhost:8080/api/logs/2601/Sample_Test
```

## 문제 해결

### 컨테이너가 시작되지 않는 경우
```bash
# 컨테이너 로그 확인
docker logs dark-gene-daemon

# 컨테이너 상태 확인
docker ps -a

# 환경 변수 확인
docker exec dark-gene-daemon env
```

### 파일 시스템 접근 권한 문제
```bash
# Docker 볼륨 마운트 확인
docker inspect dark-gene-daemon | grep -A 10 Mounts

# 파일 권한 확인
docker exec dark-gene-daemon ls -la /dark_gene_pipeline/
```

### API 호출이 실패하는 경우
```bash
# 네트워크 연결 확인
docker exec dark-gene-daemon curl http://localhost:8080/health

# Portal API 연결 테스트
docker exec dark-gene-daemon curl -v https://your-portal-domain.com/api
```

## 업데이트 방법

### 코드 변경 후 재배포
```bash
cd /home/ken/dark_gene_pipeline/daemon

# 기존 컨테이너 중지 및 삭제
docker-compose down

# 이미지 재빌드
docker-compose build --no-cache

# 재시작
docker-compose up -d

# 로그 확인
docker-compose logs -f
```

## Production 배포 시 고려사항

### 1. 보안
- `.env` 파일은 절대 Git에 커밋하지 마세요
- API Key는 주기적으로 갱신하세요
- 가능하면 HTTPS로 Portal과 통신하세요

### 2. 모니터링
- Docker 컨테이너 헬스체크 설정
- 로그 파일 크기 제한 설정
- 알림 시스템 구축 (실패 시 알림)

### 3. 백업
- 분석 결과 파일 정기 백업
- 로그 파일 보관 정책 수립

### 4. 성능
- 동시 분석 샘플 수에 따라 리소스 할당 조정
- 디스크 공간 모니터링
- 네트워크 대역폭 고려

## 시스템 서비스로 등록 (옵션)

Docker Compose를 시스템 부팅 시 자동 시작하려면:

```bash
# systemd 서비스 파일 생성
sudo nano /etc/systemd/system/dark-gene-daemon.service
```

```ini
[Unit]
Description=Dark Gene Pipeline Daemon
Requires=docker.service
After=docker.service

[Service]
Type=oneshot
RemainAfterExit=yes
WorkingDirectory=/home/ken/dark_gene_pipeline/daemon
ExecStart=/usr/bin/docker-compose up -d
ExecStop=/usr/bin/docker-compose down
TimeoutStartSec=0

[Install]
WantedBy=multi-user.target
```

```bash
# 서비스 활성화 및 시작
sudo systemctl enable dark-gene-daemon.service
sudo systemctl start dark-gene-daemon.service

# 상태 확인
sudo systemctl status dark-gene-daemon.service
```
