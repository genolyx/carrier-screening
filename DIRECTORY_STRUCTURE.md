# Dark Gene Pipeline - Directory Structure Guide

## 📁 새로운 디렉토리 구조

```
dark_gene_pipeline/
├── fastq/<work_dir>/<sample_name>/          # FASTQ 원본 파일 저장소
│   ├── *_R1_*.fq.gz                         # Read 1
│   ├── *_R2_*.fq.gz                         # Read 2
│   └── analysis.completed                   # 분석 완료 마커 (자동 생성)
│
├── analysis/<work_dir>/<sample_name>/       # 중간 파일 (intermediate files)
│   ├── pipeline_info/
│   │   ├── trace.txt                        # Task 실행 로그
│   │   ├── timeline.html                    # 타임라인 리포트
│   │   └── report.html                      # 실행 리포트
│   ├── alignment/                           # BAM 파일
│   ├── variant/                             # VCF 파일
│   ├── cnv/                                 # CNV 결과
│   ├── sv/                                  # SV 결과
│   ├── pseudogene/                          # 유사유전자 분석
│   ├── repeat/                              # 반복서열 확장
│   └── ...                                  # 기타 중간 결과
│
├── output/<work_dir>/<sample_name>/         # Portal용 최종 결과 (자동 복사)
│   ├── summary/                             # 요약 리포트
│   │   ├── *_summary_report.txt
│   │   └── *_detailed_report.txt
│   ├── snapshots/                           # 시각화
│   │   └── *_visual_report.html
│   ├── vcf/                                 # Filtered VCF
│   ├── bam/                                 # Final BAM
│   ├── cnv/                                 # CNV VCF
│   ├── sv/                                  # SV VCF
│   ├── repeat/                              # Repeat expansion
│   └── pseudogene/                          # Pseudogene results
│
└── log/<work_dir>/<sample_name>/            # 파이프라인 실행 로그
    └── nextflow.log                         # Nextflow 메인 로그
```

## 🔑 핵심 개념

### 1. Work Directory (YYMM 형식)
- 현재 날짜 기준으로 자동 생성 (예: 2601 = 2026년 1월)
- 월별로 샘플 관리 가능
- Dashboard에서 수동 생성 가능

### 2. Sample Name
- 각 샘플은 고유한 디렉토리를 가짐
- R1/R2 FASTQ 파일을 포함
- 분석 완료 시 `analysis.completed` 마커 자동 생성

### 3. 분석 완료 마커 (`analysis.completed`)
- 분석 성공 시 자동으로 fastq 디렉토리에 생성
- 이 파일이 있으면 Dashboard에서 기본적으로 숨김 처리
- "전체보기" 토글로 다시 표시 가능
- 재분석 시 자동으로 삭제됨

## 🚀 사용 방법

### 1. Dashboard 접속
```bash
cd dashboard
python3 app.py
# 브라우저: http://localhost:5000
```

### 2. 샘플 선택 및 분석
1. **Samples** 탭에서 샘플 디렉토리 선택
2. R1/R2 페어 자동 검증
3. "Analyze" 또는 "Start Selected" 클릭
4. **Monitor** 탭에서 실시간 진행 상황 확인

### 3. 결과 확인
- **Results** 탭에서 리포트 다운로드
- `output/<work_dir>/<sample_name>/` 에서 직접 접근 가능

## 📤 FASTQ 업로드

### Dashboard 업로드
1. **Upload** 탭 선택
2. Work Directory 입력 (예: 2601)
3. Sample Name 입력 (예: Sample_A10)
4. R1/R2 파일 선택 및 업로드

### 수동 업로드
```bash
# 디렉토리 생성
mkdir -p fastq/2601/Sample_A10

# FASTQ 파일 복사
cp /path/to/*_R1_*.fq.gz fastq/2601/Sample_A10/
cp /path/to/*_R2_*.fq.gz fastq/2601/Sample_A10/
```

## 🔄 재분석

완료된 샘플을 재분석하려면:

1. **전체보기** 토글 활성화
2. 완료된 샘플 선택
3. "Analyze" 클릭
4. `analysis.completed` 마커 자동 삭제 및 재분석 시작

## 🧹 정리 (Cleanup)

### 자동 정리 (권장)
파이프라인이 자동으로 처리:
- `analysis/` 중간 파일 유지 (재분석 가능)
- `output/` 최종 결과 자동 복사
- `work/` Nextflow 작업 디렉토리 삭제 (--cleanup 플래그 시)

### 수동 정리
```bash
# 특정 샘플의 중간 파일 삭제
rm -rf analysis/2601/Sample_A10/

# 전체 분석 디렉토리 정리 (주의!)
rm -rf analysis/2601/
```

## 📊 디렉토리별 용도

| 디렉토리 | 용도 | 보관 기간 | 크기 |
|---------|------|----------|------|
| **fastq/** | 원본 데이터 | 영구 | ~6GB/샘플 |
| **analysis/** | 중간 파일 | 필요시 삭제 가능 | ~15GB/샘플 |
| **output/** | 최종 결과 | 영구 | ~2GB/샘플 |
| **log/** | 실행 로그 | 1개월 | ~10MB/샘플 |

## 🔍 문제 해결

### Q: 샘플이 목록에 나타나지 않아요
A: 
1. `fastq/<work_dir>/<sample_name>/` 구조가 올바른지 확인
2. R1/R2 파일명에 `_R1`, `_R2` 포함 확인
3. "전체보기" 토글 활성화 (완료된 샘플인 경우)

### Q: "Missing Pair" 에러가 발생해요
A:
- R1과 R2 파일이 같은 디렉토리에 있는지 확인
- 파일명에 `_R1_`, `_R2_` 또는 `_1.fq`, `_2.fq` 포함 확인

### Q: 재분석이 안 돼요
A:
- "전체보기"로 완료된 샘플 표시
- 선택 후 "Analyze" 클릭 시 `analysis.completed` 자동 삭제

## 🔐 권한 설정

```bash
# 디렉토리 권한 확인
chmod 755 fastq/ analysis/ output/ log/

# 샘플 디렉토리 권한
chmod 755 fastq/2601/Sample_*
```

## 📝 마이그레이션

기존 `fastq_old/` 에서 새 구조로 마이그레이션:

```bash
./migrate_fastq.sh
```

스크립트가 자동으로:
1. 샘플명 추출
2. 디렉토리 생성
3. 파일 복사
4. `analysis.completed` 마커 생성

## 🎯 Best Practices

1. **Work Directory 명명 규칙**: YYMM 형식 권장 (2601, 2602, ...)
2. **Sample Name**: 명확하고 고유한 이름 사용
3. **정기 정리**: 월 1회 `analysis/` 디렉토리 정리
4. **백업**: `fastq/`와 `output/` 디렉토리는 정기적으로 백업

## 📞 지원

문제가 발생하면 로그 확인:
```bash
# Nextflow 로그
cat log/<work_dir>/<sample_name>/nextflow.log

# Pipeline trace
cat analysis/<work_dir>/<sample_name>/pipeline_info/trace.txt
```
