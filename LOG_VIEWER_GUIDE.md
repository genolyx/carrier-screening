# Nextflow 로그 확인 가이드

## 🔍 로그 확인 방법

### 1. 서버에서 직접 확인

#### A. 스크립트 사용 (권장)

```bash
# 로그 목록 보기
./scripts/view_logs.sh --list

# 특정 샘플 로그 보기
./scripts/view_logs.sh -w 2601 -s Sample_A10

# 마지막 100줄만 보기
./scripts/view_logs.sh -w 2601 -s Sample_A10 -n 100

# 실시간 모니터링 (tail -f)
./scripts/view_logs.sh -w 2601 -s Sample_A10 -f

# 최신 로그 자동 선택하여 실시간 모니터링
./scripts/view_logs.sh --latest -f
```

#### B. 직접 파일 접근

```bash
# 로그 파일 위치
log/<work_dir>/<sample_name>/nextflow.log

# 예시
cat log/2601/Sample_A10/nextflow.log

# tail -f로 실시간 모니터링
tail -f log/2601/Sample_A10/nextflow.log

# 마지막 100줄
tail -n 100 log/2601/Sample_A10/nextflow.log

# 에러만 필터링
grep -i error log/2601/Sample_A10/nextflow.log

# 특정 프로세스 로그
grep "ALIGN_AND_SORT" log/2601/Sample_A10/nextflow.log
```

---

### 2. Portal에서 확인

#### A. Daemon API를 통한 로그 조회

Portal 백엔드에서 호출:

```python
import requests

# 1. 로그 목록 조회
response = requests.get('http://daemon-server:8080/api/logs/list')
logs = response.json()['logs']

# 2. 특정 샘플 로그 전체 조회
response = requests.get('http://daemon-server:8080/api/logs/2601/Sample_A10')
log_content = response.json()['content']

# 3. 마지막 100줄만 조회
response = requests.get('http://daemon-server:8080/api/logs/2601/Sample_A10?tail=100')
log_content = response.json()['content']

# 4. 로그 파일 다운로드
response = requests.get('http://daemon-server:8080/api/logs/2601/Sample_A10/download')
with open('nextflow.log', 'wb') as f:
    f.write(response.content)
```

#### B. 실시간 로그 스트리밍 (SSE)

Portal 프론트엔드에서:

```javascript
// EventSource로 실시간 로그 수신
const eventSource = new EventSource(
    'http://daemon-server:8080/api/logs/2601/Sample_A10?follow=true&tail=100'
);

eventSource.onmessage = function(event) {
    const logLine = event.data;
    
    // 로그 라인을 화면에 추가
    const logContainer = document.getElementById('log-container');
    logContainer.innerHTML += logLine + '\n';
    
    // Auto-scroll
    logContainer.scrollTop = logContainer.scrollHeight;
};

eventSource.onerror = function(err) {
    console.error('EventSource error:', err);
    eventSource.close();
};

// 중지하려면
eventSource.close();
```

#### C. HTML Log Viewer (독립 뷰어)

```bash
# 브라우저에서 열기
open daemon/log_viewer.html?work_dir=2601&sample_name=Sample_A10&daemon_url=http://localhost:8080
```

**기능:**
- ✅ 실시간 tail -f 모드
- ✅ 구문 강조 (ERROR, WARNING, SUCCESS)
- ✅ 검색 기능
- ✅ Auto-scroll
- ✅ 다운로드
- ✅ 라인 수 선택 (100, 500, 1000, All)

---

## 📡 Daemon API 엔드포인트

### 1. 로그 목록 조회
```http
GET /api/logs/list

Response:
{
    "logs": [
        {
            "work_dir": "2601",
            "sample_name": "Sample_A10",
            "log_file": "/pipeline/log/2601/Sample_A10/nextflow.log",
            "size": 12345,
            "modified": "2026-01-21T14:30:00"
        }
    ],
    "total": 1
}
```

### 2. 로그 내용 조회
```http
GET /api/logs/{work_dir}/{sample_name}?tail=100

Response:
{
    "work_dir": "2601",
    "sample_name": "Sample_A10",
    "log_file": "/pipeline/log/2601/Sample_A10/nextflow.log",
    "content": "...",
    "lines": 100
}
```

### 3. 실시간 로그 스트리밍
```http
GET /api/logs/{work_dir}/{sample_name}?follow=true&tail=100

Response: (Server-Sent Events)
data: [2026-01-21 14:30:00] INFO  Starting analysis...
data: [2026-01-21 14:30:01] INFO  Aligning reads...
...
```

### 4. 로그 파일 다운로드
```http
GET /api/logs/{work_dir}/{sample_name}/download

Response: (File download)
Content-Disposition: attachment; filename=2601_Sample_A10_nextflow.log
```

---

## 🌐 Portal 통합 예시

### Portal Summary 페이지에 로그 버튼 추가

```html
<!-- Portal: /analysis/summary 페이지 -->
<table>
    <thead>
        <tr>
            <th>Order ID</th>
            <th>Sample</th>
            <th>Status</th>
            <th>Progress</th>
            <th>Actions</th>
        </tr>
    </thead>
    <tbody>
        <tr>
            <td>ORD-2601-001</td>
            <td>2601/Sample_A10</td>
            <td><span class="badge badge-running">RUNNING</span></td>
            <td>45%</td>
            <td>
                <button onclick="viewLogs('2601', 'Sample_A10')">
                    📋 View Logs
                </button>
            </td>
        </tr>
    </tbody>
</table>

<script>
function viewLogs(workDir, sampleName) {
    // 모달 또는 새 창에서 로그 표시
    const url = `/analysis/logs?work_dir=${workDir}&sample_name=${sampleName}`;
    window.open(url, '_blank', 'width=1200,height=800');
}
</script>
```

### Portal 로그 페이지 구현

```python
# Portal Backend (Flask/Django)
from flask import Flask, render_template, request
import requests

app = Flask(__name__)
DAEMON_URL = 'http://daemon-server:8080'

@app.route('/analysis/logs')
def view_logs():
    work_dir = request.args.get('work_dir')
    sample_name = request.args.get('sample_name')
    
    # Daemon에서 로그 조회
    response = requests.get(f'{DAEMON_URL}/api/logs/{work_dir}/{sample_name}?tail=500')
    log_data = response.json()
    
    return render_template('logs.html', 
        work_dir=work_dir,
        sample_name=sample_name,
        log_content=log_data['content'],
        daemon_url=DAEMON_URL
    )
```

```html
<!-- Portal: templates/logs.html -->
<!DOCTYPE html>
<html>
<head>
    <title>Nextflow Logs - {{ work_dir }}/{{ sample_name }}</title>
    <style>
        /* 위의 log_viewer.html 스타일 재사용 */
    </style>
</head>
<body>
    <div class="container">
        <h1>Nextflow Logs: {{ work_dir }}/{{ sample_name }}</h1>
        
        <div class="controls">
            <button id="btn-follow" onclick="toggleFollow()">▶ Start Following</button>
            <button onclick="downloadLog()">⬇️ Download</button>
        </div>
        
        <pre id="log-container">{{ log_content }}</pre>
    </div>
    
    <script>
        const workDir = '{{ work_dir }}';
        const sampleName = '{{ sample_name }}';
        const daemonUrl = '{{ daemon_url }}';
        
        let eventSource = null;
        let isFollowing = false;
        
        function toggleFollow() {
            if (isFollowing) {
                eventSource.close();
                isFollowing = false;
                document.getElementById('btn-follow').textContent = '▶ Start Following';
            } else {
                const url = `${daemonUrl}/api/logs/${workDir}/${sampleName}?follow=true&tail=100`;
                eventSource = new EventSource(url);
                
                eventSource.onmessage = function(event) {
                    const container = document.getElementById('log-container');
                    container.textContent += event.data + '\n';
                    container.scrollTop = container.scrollHeight;
                };
                
                isFollowing = true;
                document.getElementById('btn-follow').textContent = '⏸ Stop Following';
            }
        }
        
        function downloadLog() {
            window.location.href = `${daemonUrl}/api/logs/${workDir}/${sampleName}/download`;
        }
    </script>
</body>
</html>
```

---

## 🔧 로그 분석 유틸리티

### 에러 요약
```bash
# 에러만 추출
grep -i "error\|failed" log/2601/Sample_A10/nextflow.log

# 에러 카운트
grep -i "error" log/2601/Sample_A10/nextflow.log | wc -l
```

### 진행 상황 파악
```bash
# 완료된 프로세스
grep "Completed" log/2601/Sample_A10/nextflow.log

# 현재 실행 중인 프로세스
grep "Submitted\|Running" log/2601/Sample_A10/nextflow.log | tail -10

# Duration 확인
grep "Duration" log/2601/Sample_A10/nextflow.log
```

### 리소스 사용량
```bash
# CPU/Memory 정보
grep "%cpu\|%mem" log/2601/Sample_A10/nextflow.log

# Peak memory
grep "peak_rss\|peak_vmem" log/2601/Sample_A10/nextflow.log
```

---

## 📊 로그 포맷

Nextflow 로그는 다음 형식으로 출력됩니다:

```
Jan-21 14:30:00.123 [main] INFO  nextflow.cli.Launcher - Starting Nextflow...
Jan-21 14:30:01.456 [Task monitor] DEBUG nextflow.processor.TaskProcessor - [warm up] executor > local
Jan-21 14:30:05.789 [Task submitter] INFO  nextflow.Session - [a1/b2c3d4] Submitted process > ALIGN_AND_SORT (Sample_A10)
Jan-21 14:35:10.123 [Task monitor] INFO  nextflow.Session - [a1/b2c3d4] Completed process > ALIGN_AND_SORT (Sample_A10)
...
Jan-21 16:45:00.000 [main] INFO  nextflow.cli.CmdRun - Execution complete -- Goodbye
```

**주요 키워드:**
- `Submitted`: 프로세스 제출
- `Completed`: 프로세스 완료
- `ERROR`: 에러 발생
- `WARNING`: 경고
- `Duration`: 총 실행 시간
- `Execution complete -- Goodbye`: 성공적 완료

---

## 🎯 Portal 구현 체크리스트

- [ ] Daemon API 엔드포인트 테스트
- [ ] Portal에 "View Logs" 버튼 추가
- [ ] 로그 페이지 구현
- [ ] 실시간 스트리밍 구현 (optional)
- [ ] 로그 다운로드 기능
- [ ] 에러 하이라이팅
- [ ] 검색 기능
- [ ] Auto-scroll

---

## 🚀 빠른 시작

### 1. Daemon에 로그 API 추가 (완료)
```bash
# Daemon 재시작
cd daemon
docker-compose restart
```

### 2. 서버에서 로그 확인
```bash
# 스크립트 사용
./scripts/view_logs.sh --list
./scripts/view_logs.sh --latest -f
```

### 3. Portal에서 로그 확인
```bash
# API 테스트
curl http://localhost:8080/api/logs/list

# 특정 로그 조회
curl http://localhost:8080/api/logs/2601/Sample_A10?tail=100
```

### 4. HTML 뷰어 사용
```bash
# 브라우저에서
open daemon/log_viewer.html?work_dir=2601&sample_name=Sample_A10
```

모든 준비가 완료되었습니다! 🎉
