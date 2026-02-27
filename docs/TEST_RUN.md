# 테스트 실행 가이드

## 1. Dashboard 시작

```bash
cd /home/ken/dark_gene_pipeline/dashboard
python3 app.py
```

브라우저에서 http://localhost:5000 접속

## 2. 현재 상태 확인

### 샘플 목록
- `2601/Exome2Comp_NA12878_LPEF2_TwistCap_G10_S107` (완료됨, R1/R2 존재)
- `2601/Sample_Test` (미완료, R1/R2 없음)

### 디렉토리 구조
```
fastq/2601/
├── Exome2Comp_NA12878_LPEF2_TwistCap_G10_S107/
│   ├── analysis.completed
│   ├── *_R1_*.fq.gz
│   └── *_R2_*.fq.gz
└── Sample_Test/ (빈 디렉토리)
```

## 3. 테스트 시나리오

### 시나리오 1: 완료된 샘플 숨김/표시
1. Dashboard 접속
2. 기본 상태: Sample_Test만 표시 (완료된 샘플 숨김)
3. "Show Completed" 토글 활성화
4. 확인: 2개 샘플 모두 표시

### 시나리오 2: 완료된 샘플 재분석
1. "Show Completed" 토글 ON
2. "Exome2Comp..." 샘플 선택
3. "Analyze" 클릭
4. 확인: analysis.completed 삭제됨
5. Monitor 탭에서 진행 상황 확인

### 시나리오 3: 새 샘플 업로드
1. Upload 탭 선택
2. Work Directory: 2601
3. Sample Name: Sample_New
4. R1/R2 파일 선택 후 업로드
5. Samples 탭으로 돌아가서 확인

## 4. 파이프라인 실행 테스트 (Dry Run)

```bash
cd /home/ken/dark_gene_pipeline

# 테스트 실행 (실제 분석 없이 구조만 확인)
./nextflow run main.nf \
  --fastq_dir fastq/2601/Exome2Comp_NA12878_LPEF2_TwistCap_G10_S107 \
  --outdir analysis/2601/Exome2Comp_NA12878_LPEF2_TwistCap_G10_S107 \
  --output_dir output/2601/Exome2Comp_NA12878_LPEF2_TwistCap_G10_S107 \
  --sample_name Exome2Comp_NA12878_LPEF2_TwistCap_G10_S107 \
  -preview
```

## 5. 예상 결과

### 분석 완료 후 디렉토리 구조
```
analysis/2601/Exome2Comp.../
├── pipeline_info/
│   ├── trace.txt
│   ├── timeline.html
│   └── report.html
├── alignment/
│   └── *.md.bam
├── variant/
│   └── *_filtered.vcf.gz
└── ...

output/2601/Exome2Comp.../
├── summary/
│   ├── *_summary_report.txt
│   └── *_detailed_report.txt
├── snapshots/
│   └── *_visual_report.html
├── vcf/
├── bam/
└── ...

log/2601/Exome2Comp.../
└── nextflow.log
```

## 6. 검증 체크리스트

- [ ] Dashboard 정상 시작
- [ ] 샘플 목록 표시 (미완료만)
- [ ] "전체보기" 토글 작동
- [ ] 완료된 샘플 재선택 가능
- [ ] R1/R2 검증 작동
- [ ] 분석 시작 가능
- [ ] Monitor 탭 진행률 표시
- [ ] 결과 파일 output 디렉토리에 복사
- [ ] analysis.completed 자동 생성

## 7. 로그 확인

```bash
# Dashboard 로그
tail -f /home/ken/dark_gene_pipeline/dashboard.log

# Nextflow 로그
tail -f log/2601/*/nextflow.log

# Trace 확인
cat analysis/2601/*/pipeline_info/trace.txt
```
