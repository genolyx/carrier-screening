#!/usr/bin/env python3
"""
Daemon API Server

Portal에서 호출할 수 있는 REST API 제공:
- GET  /api/status            : 전체 상태 조회
- GET  /api/summary           : 분석 요약 (첨부 이미지와 같은 형식)
- GET  /api/orders            : 주문 목록
- POST /api/orders/{id}/retry : 재시도
"""

from flask import Flask, jsonify, request
from daemon import AnalysisMonitor, PortalAPI
import os

app = Flask(__name__)

# Configuration
PORTAL_URL = os.environ.get('PORTAL_URL', 'https://portal.genolyx.com')
API_KEY = os.environ.get('PORTAL_API_KEY', '')

portal_api = PortalAPI(PORTAL_URL, API_KEY)
monitor = AnalysisMonitor(portal_api)


@app.route('/api/health', methods=['GET'])
def health():
    """헬스 체크"""
    return jsonify({'status': 'healthy', 'service': 'dark-gene-daemon'})


@app.route('/api/summary', methods=['GET'])
def get_summary():
    """
    분석 상태 요약 (Portal Summary 페이지용)
    
    Response:
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
    """
    summary = monitor.get_summary()
    return jsonify(summary)


@app.route('/api/status', methods=['GET'])
def get_status():
    """
    전체 상태 조회
    
    Response:
    {
        "daemon": "running",
        "samples": {...},
        "summary": {...}
    }
    """
    samples = monitor.scan_samples()
    summary = monitor.get_summary()
    
    return jsonify({
        'daemon': 'running',
        'samples': samples,
        'summary': summary,
        'timestamp': monitor.state.get('last_updated', '')
    })


@app.route('/api/orders', methods=['GET'])
def get_orders():
    """
    주문 목록 조회
    
    Query Params:
    - status: RUNNING, COMPLETED, FAILED
    - date: YYYY-MM-DD
    """
    status_filter = request.args.get('status')
    date_filter = request.args.get('date')
    
    orders = []
    samples = monitor.scan_samples()
    
    for sample in samples:
        sample_key = f"{sample['work_dir']}/{sample['sample_name']}"
        order_info = monitor.state['orders'].get(sample_key, {})
        
        if not order_info:
            continue
        
        # 상태 필터
        if status_filter:
            if sample['is_running'] and status_filter != 'RUNNING':
                continue
            if sample['is_completed'] and status_filter != 'COMPLETED':
                continue
        
        # 날짜 필터
        if date_filter:
            created_at = order_info.get('created_at', '')
            if not created_at.startswith(date_filter):
                continue
        
        # 진행률
        progress = 100 if sample['is_completed'] else monitor.get_analysis_progress(
            sample['work_dir'], sample['sample_name']
        )
        
        orders.append({
            'order_id': order_info.get('order_id'),
            'work_dir': sample['work_dir'],
            'sample_name': sample['sample_name'],
            'status': 'COMPLETED' if sample['is_completed'] else 'RUNNING' if sample['is_running'] else 'WAITING',
            'progress': progress,
            'created_at': order_info.get('created_at'),
            'uploaded_at': order_info.get('uploaded_at'),
            'uploaded': order_info.get('uploaded', False)
        })
    
    return jsonify({'orders': orders, 'total': len(orders)})


@app.route('/api/orders/submit', methods=['POST'])
def submit_order():
    """
    Portal에서 Order 제출 (Use Case 2)
    
    Request Body:
    {
        "order_id": "ORD-2026-001",
        "work_dir": "2601",
        "sample_name": "Sample_A10",
        "fastq_r1_url": "https://portal.com/files/xxx_R1.fastq.gz",
        "fastq_r2_url": "https://portal.com/files/xxx_R2.fastq.gz",
        "priority": "normal",  // "urgent", "normal", "low"
        "notify_email": "user@example.com"
    }
    
    Response:
    {
        "status": "accepted",
        "order_id": "ORD-2026-001",
        "message": "Order submitted successfully"
    }
    """
    import os
    import subprocess
    from datetime import datetime
    
    data = request.json
    
    # 필수 파라미터 검증
    required_fields = ['order_id', 'work_dir', 'sample_name', 'fastq_r1_url', 'fastq_r2_url']
    for field in required_fields:
        if field not in data:
            return jsonify({'error': f'Missing required field: {field}'}), 400
    
    order_id = data['order_id']
    work_dir = data['work_dir']
    sample_name = data['sample_name']
    r1_url = data['fastq_r1_url']
    r2_url = data['fastq_r2_url']
    
    # FASTQ 디렉토리 생성
    BASE_DIR = os.environ.get('BASE_DIR', '/data')
    fastq_dir = os.path.join(BASE_DIR, 'fastq', work_dir, sample_name)
    os.makedirs(fastq_dir, exist_ok=True)
    
    # R1 파일명 추출
    r1_filename = os.path.basename(r1_url)
    r2_filename = os.path.basename(r2_url)
    
    r1_path = os.path.join(fastq_dir, r1_filename)
    r2_path = os.path.join(fastq_dir, r2_filename)
    
    # FASTQ 파일 다운로드
    try:
        # wget 또는 curl로 다운로드
        subprocess.run(['wget', '-O', r1_path, r1_url], check=True, capture_output=True)
        subprocess.run(['wget', '-O', r2_path, r2_url], check=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        return jsonify({
            'error': 'Failed to download FASTQ files',
            'details': str(e)
        }), 500
    
    # Order 정보 저장
    sample_key = f"{work_dir}/{sample_name}"
    monitor.state['orders'][sample_key] = {
        'order_id': order_id,
        'created_at': datetime.now().isoformat(),
        'status': 'WAITING',
        'progress': 0,
        'uploaded': False,
        'priority': data.get('priority', 'normal'),
        'notify_email': data.get('notify_email', '')
    }
    monitor.save_state()
    
    # Portal에 접수 확인 전송
    monitor.portal_api.update_order_status(order_id, 'WAITING', 0, 'Order accepted')
    
    # Dashboard API 호출하여 분석 시작 (또는 큐에 추가)
    # 여기서는 자동으로 분석을 시작하지 않고, Dashboard에서 수동 시작하도록 함
    # 자동 시작을 원하면 아래 주석 해제:
    # try:
    #     import requests
    #     requests.post('http://localhost:5000/start', json={
    #         'work_dir': work_dir,
    #         'sample_name': sample_name
    #     })
    # except Exception as e:
    #     pass  # 실패해도 Order는 접수됨
    
    return jsonify({
        'status': 'accepted',
        'order_id': order_id,
        'work_dir': work_dir,
        'sample_name': sample_name,
        'message': 'Order submitted successfully. Waiting for analysis to start.'
    }), 202


@app.route('/api/orders/<order_id>/retry', methods=['POST'])
def retry_order(order_id):
    """주문 재시도"""
    # Order ID로 샘플 찾기
    for sample_key, order_info in monitor.state['orders'].items():
        if order_info.get('order_id') == order_id:
            work_dir, sample_name = sample_key.split('/')
            
            # analysis.completed 삭제
            completed_marker = os.path.join(
                monitor.state.get('FASTQ_BASE', '/pipeline/fastq'),
                work_dir,
                sample_name,
                'analysis.completed'
            )
            
            if os.path.exists(completed_marker):
                os.remove(completed_marker)
            
            # 상태 초기화
            order_info['uploaded'] = False
            order_info['uploaded_at'] = None
            monitor.save_state()
            
            # Portal 상태 업데이트
            monitor.portal_api.update_analysis_status(order_id, 'WAITING', 0)
            
            return jsonify({'message': 'Order reset for retry', 'order_id': order_id})
    
    return jsonify({'error': 'Order not found'}), 404


@app.route('/api/samples/scan', methods=['POST'])
def scan_samples():
    """샘플 재스캔 (수동 트리거)"""
    monitor.process_completed_samples()
    monitor.process_running_samples()
    return jsonify({'message': 'Scan completed'})


@app.route('/api/logs/<work_dir>/<sample_name>', methods=['GET'])
def get_logs(work_dir, sample_name):
    """
    Nextflow 로그 조회
    
    Query Params:
    - tail: 마지막 N줄만 (default: 전체)
    - follow: true일 경우 실시간 스트리밍 (SSE)
    
    Example:
    - GET /api/logs/2601/Sample_A10?tail=100
    - GET /api/logs/2601/Sample_A10?follow=true
    """
    import os
    from flask import Response, stream_with_context
    
    LOG_BASE = os.environ.get('PIPELINE_BASE_DIR', '/pipeline') + '/log'
    log_file = os.path.join(LOG_BASE, work_dir, sample_name, 'nextflow.log')
    
    if not os.path.exists(log_file):
        return jsonify({'error': 'Log file not found', 'path': log_file}), 404
    
    # Query params
    tail = request.args.get('tail', type=int)
    follow = request.args.get('follow', 'false').lower() == 'true'
    
    if follow:
        # SSE (Server-Sent Events) 스트리밍
        def generate():
            import time
            import subprocess
            
            # tail -f 사용
            proc = subprocess.Popen(
                ['tail', '-f', '-n', str(tail or 100), log_file],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
            
            try:
                for line in iter(proc.stdout.readline, ''):
                    if line:
                        yield f"data: {line}\n\n"
                    time.sleep(0.1)
            finally:
                proc.kill()
        
        return Response(
            stream_with_context(generate()),
            mimetype='text/event-stream',
            headers={
                'Cache-Control': 'no-cache',
                'X-Accel-Buffering': 'no'
            }
        )
    
    else:
        # 일반 로그 반환
        try:
            with open(log_file, 'r') as f:
                if tail:
                    # 마지막 N줄만
                    lines = f.readlines()
                    content = ''.join(lines[-tail:])
                else:
                    # 전체
                    content = f.read()
            
            return jsonify({
                'work_dir': work_dir,
                'sample_name': sample_name,
                'log_file': log_file,
                'content': content,
                'lines': len(content.splitlines())
            })
        except Exception as e:
            return jsonify({'error': str(e)}), 500


@app.route('/api/logs/<work_dir>/<sample_name>/download', methods=['GET'])
def download_log(work_dir, sample_name):
    """로그 파일 다운로드"""
    import os
    from flask import send_file
    
    LOG_BASE = os.environ.get('PIPELINE_BASE_DIR', '/pipeline') + '/log'
    log_file = os.path.join(LOG_BASE, work_dir, sample_name, 'nextflow.log')
    
    if not os.path.exists(log_file):
        return jsonify({'error': 'Log file not found'}), 404
    
    return send_file(
        log_file,
        as_attachment=True,
        download_name=f'{work_dir}_{sample_name}_nextflow.log'
    )


@app.route('/api/logs/list', methods=['GET'])
def list_logs():
    """
    모든 로그 파일 목록
    
    Response:
    [
        {
            "work_dir": "2601",
            "sample_name": "Sample_A10",
            "log_file": "/pipeline/log/2601/Sample_A10/nextflow.log",
            "size": 12345,
            "modified": "2026-01-21T14:30:00"
        }
    ]
    """
    import os
    from datetime import datetime
    
    LOG_BASE = os.environ.get('PIPELINE_BASE_DIR', '/pipeline') + '/log'
    logs = []
    
    if not os.path.exists(LOG_BASE):
        return jsonify({'logs': []})
    
    for work_dir in os.listdir(LOG_BASE):
        work_path = os.path.join(LOG_BASE, work_dir)
        if not os.path.isdir(work_path):
            continue
        
        for sample_name in os.listdir(work_path):
            sample_path = os.path.join(work_path, sample_name)
            if not os.path.isdir(sample_path):
                continue
            
            log_file = os.path.join(sample_path, 'nextflow.log')
            if os.path.exists(log_file):
                stat = os.stat(log_file)
                logs.append({
                    'work_dir': work_dir,
                    'sample_name': sample_name,
                    'log_file': log_file,
                    'size': stat.st_size,
                    'modified': datetime.fromtimestamp(stat.st_mtime).isoformat()
                })
    
    return jsonify({'logs': logs, 'total': len(logs)})


if __name__ == '__main__':
    app.run(host='0.0.0.0', port=8080, debug=False)
