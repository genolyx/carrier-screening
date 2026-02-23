# 최종 디렉토리 구조 완성 요약

## ✅ 완료된 변경 사항

### 1. publishDir 모드 통일
모든 프로세스를 `mode: 'copy'`로 변경하여 work/ 삭제 후에도 안전하게 보관

```diff
# modules/align.nf
- publishDir "${params.outdir}/alignment", mode: 'symlink'
+ publishDir "${params.outdir}/alignment", mode: 'copy'

# modules/cnv.nf (GCNV_CLARITY)
- publishDir "${params.outdir}/cnv", mode: 'symlink'
+ publishDir "${params.outdir}/cnv", mode: 'copy'
```

### 2. 디렉토리 구조 확정

```
dark_gene_pipeline/
│
├── fastq/<work_dir>/<sample_name>/
│   ├── *_R1_*.fq.gz                    # 원본 FASTQ
│   ├── *_R2_*.fq.gz
│   └── analysis.completed              # 완료 마커
│
├── analysis/<work_dir>/<sample_name>/  # 모든 intermediate 파일
│   ├── alignment/
│   │   ├── *.bam                       # BAM 파일 (copy)
│   │   └── *.bam.bai
│   ├── variant/
│   │   └── *_filtered.vcf.gz           # VCF 파일 (copy)
│   ├── cnv/
│   │   ├── counts/                     # Read counts
│   │   ├── segments/                   # CNV segments
│   │   └── cohort/gcnv-model/          # 재사용 가능한 모델
│   ├── sv/
│   │   └── *_manta.vcf.gz              # SV calls
│   ├── pseudogene/
│   │   ├── paraphase/                  # Paraphase 결과
│   │   └── smaca/                      # SMAca 결과
│   ├── repeat/
│   │   ├── *_eh.vcf                    # ExpansionHunter
│   │   ├── *_eh.json
│   │   └── *.svg                       # REViewer 시각화
│   ├── fallback/
│   │   ├── *_hba_fallback.txt
│   │   └── *_cyp21a2_fallback.txt
│   ├── coverage/
│   │   └── *_intron_depth.txt
│   ├── snapshots/
│   │   └── *_visual_report.html
│   ├── summary/
│   │   ├── *_summary_report.txt
│   │   └── *_detailed_report.txt
│   └── pipeline_info/
│       ├── trace.txt
│       ├── timeline.html
│       └── report.html
│
├── output/<work_dir>/<sample_name>/    # 포털용 최종 결과만
│   ├── summary/
│   │   ├── *_summary_report.txt        # 복사
│   │   └── *_detailed_report.txt       # 복사
│   ├── snapshots/
│   │   └── *_visual_report.html        # 복사
│   ├── vcf/
│   │   └── *_filtered.vcf.gz           # 복사
│   ├── bam/
│   │   └── *.md.bam                    # 복사
│   ├── cnv/
│   │   └── *_cnv.vcf.gz                # 복사
│   ├── sv/
│   │   └── *_manta.vcf.gz              # 복사
│   ├── repeat/
│   │   ├── *.vcf, *.json, *.svg        # 복사
│   └── pseudogene/
│       └── paraphase/, smaca/          # 복사
│
├── log/<work_dir>/<sample_name>/
│   └── nextflow.log                    # Nextflow 실행 로그
│
└── work/                                # Nextflow 작업 디렉토리
    └── [해시]/                         # 분석 완료 후 삭제됨
```

---

## 🔄 파일 흐름

### 분석 중
```
1. FASTQ 입력
   fastq/2601/Sample_A10/*_R1_*.fq.gz

2. Nextflow 작업
   work/a1/b2c3.../
   └── sample.bam (실제 파일 생성)

3. publishDir (copy)
   analysis/2601/Sample_A10/alignment/sample.bam
   (work/에서 복사됨)

4. 포털용 복사 (onComplete)
   output/2601/Sample_A10/bam/sample.bam
   (analysis/에서 복사됨)
```

### 완료 후
```
work/                           → 삭제됨 (--cleanup 시)
analysis/2601/Sample_A10/       → ✅ 보존 (재분석 가능)
output/2601/Sample_A10/         → ✅ 보존 (포털 제공)
fastq/2601/Sample_A10/          → ✅ 보존 (원본)
log/2601/Sample_A10/            → ✅ 보존 (디버깅)
```

---

## 📊 디스크 사용량

### 12 샘플 배치 기준

#### 분석 중
```
fastq/              : 72GB  (6GB/샘플 × 12)
work/               : 180GB (15GB/샘플 × 12, 임시)
analysis/           : 0GB   (아직 복사 안됨)
output/             : 0GB   (완료 후 생성)
log/                : 120MB (10MB/샘플 × 12)
-------------------------------------------
총합                : 252GB
```

#### 분석 완료 후 (--cleanup)
```
fastq/              : 72GB  (원본 보존)
work/               : 0GB   (삭제됨)
analysis/           : 60GB  (5GB/샘플 × 12, intermediate)
output/             : 24GB  (2GB/샘플 × 12, 포털용)
log/                : 120MB
-------------------------------------------
총합                : 156GB
```

#### 기존 구조 (비교)
```
fastq/              : 72GB
work/               : 0GB (삭제됨)
results/            : ⚠️ symlink 깨짐
final_archived_results/ : 40GB
-------------------------------------------
총합                : 112GB (재분석 불가)
```

### 개선 효과
- ✅ 재분석 가능 (analysis/ 보존)
- ✅ 포털 제공 (output/ 별도)
- ✅ 디스크 효율 (work/ 삭제)
- ✅ 명확한 구조 (용도별 분리)

---

## 🎯 사용 시나리오

### 시나리오 1: 일반 분석
```bash
# Dashboard에서 샘플 선택 → Analyze

# 내부 실행:
./nextflow run main.nf \
  --fastq_dir fastq/2601/Sample_A10 \
  --outdir analysis/2601/Sample_A10 \
  --output_dir output/2601/Sample_A10 \
  --sample_name Sample_A10 \
  --cleanup

# 결과:
analysis/2601/Sample_A10/  ✅ 모든 중간 파일
output/2601/Sample_A10/    ✅ 포털용 최종 결과
work/                      ❌ 삭제됨
fastq/2601/Sample_A10/analysis.completed ✅ 생성됨
```

### 시나리오 2: 재분석
```bash
# Dashboard "전체보기" → 완료된 샘플 선택 → Analyze

# 내부 동작:
1. analysis.completed 삭제
2. analysis/ 디렉토리 활용 (재생성 스킵 가능)
3. 동일한 파이프라인 실행
```

### 시나리오 3: 특정 단계만 재실행
```bash
# analysis/ 디렉토리가 있으므로 Nextflow 캐싱 활용
./nextflow run main.nf \
  --fastq_dir fastq/2601/Sample_A10 \
  --outdir analysis/2601/Sample_A10 \
  -resume

# 변경된 부분만 재실행됨
```

---

## 🔧 핵심 변경 요약

### 1. publishDir 모드
```
symlink → copy
```
→ work/ 삭제해도 analysis/ 안전

### 2. onComplete 로직
```groovy
if (workflow.success) {
    // 1. output/ 복사 (포털용)
    copyToOutputDir()
    
    // 2. work/ 삭제 (디스크 확보)
    if (params.cleanup) {
        rm -rf work/
    }
    
    // 3. analysis.completed 생성
    touch fastq/.../analysis.completed
}
```

### 3. 제거된 요소
- ❌ `final_archived_results/` (불필요, analysis/로 통합)
- ❌ `scripts/archive_results.sh` (onComplete로 대체)
- ❌ `active_run/` (fastq/ 디렉토리로 통합)
- ❌ `fastq_old/` (analysis.completed 마커로 대체)

---

## 📝 장점 요약

### 명확성
- ✅ intermediate (analysis/) vs final (output/) 명확 분리
- ✅ 용도별 디렉토리 (fastq, analysis, output, log)
- ✅ work_dir/sample_name 계층 구조

### 효율성
- ✅ work/ 삭제로 디스크 절약
- ✅ publishDir copy로 재분석 가능
- ✅ output/만 포털에 제공 (불필요한 파일 배제)

### 안정성
- ✅ work/ 삭제해도 안전 (analysis/ 보존)
- ✅ analysis.completed 마커로 중복 방지
- ✅ 재분석 시 Nextflow 캐싱 활용

### 확장성
- ✅ work_dir (YYMM) 기반 월별 관리
- ✅ 샘플별 독립 디렉토리
- ✅ 대량 샘플 동시 처리 가능

---

## 🚀 다음 단계

### 1. 테스트
```bash
cd /home/ken/dark_gene_pipeline/dashboard
python3 app.py

# 브라우저: http://localhost:5000
# TEST_RUN.md 시나리오 실행
```

### 2. 프로덕션 배포
```bash
# 기존 데이터 백업
tar -czf backup_$(date +%Y%m%d).tar.gz fastq_old/ final_archived_results/

# 마이그레이션
./migrate_fastq.sh

# Dashboard 재시작
cd dashboard && python3 app.py
```

### 3. 사용자 교육
- `DIRECTORY_STRUCTURE.md` 배포
- `WORKFLOW_ANALYSIS.md` 참고 자료 제공
- 새로운 인터페이스 교육

---

## ✅ 체크리스트

- [x] publishDir 모드 변경 (symlink → copy)
- [x] 디렉토리 구조 재설계
- [x] Dashboard 전면 수정
- [x] main.nf onComplete 로직 구현
- [x] 마이그레이션 스크립트 작성
- [x] 문서화 완료
- [ ] 실제 샘플 테스트
- [ ] 프로덕션 배포

