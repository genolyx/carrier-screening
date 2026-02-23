# Nextflow 워크플로우 분석 - 파일 생성 경로

## 📁 파일 생성 구조 분석

### 1. Nextflow의 3단계 파일 관리

```
[1단계] work/ 디렉토리 (실제 작업 공간)
   ↓
[2단계] results/ 디렉토리 (publishDir로 복사/링크)
   ↓
[3단계] final_archived_results/ (archive_results.sh로 복사)
```

---

## 🔍 상세 분석

### 1단계: work/ 디렉토리 (Nextflow 작업 공간)

#### 특징
- Nextflow가 **모든 프로세스를 실행하는 실제 공간**
- 각 task마다 고유한 해시 디렉토리 생성 (예: `work/a1/b2c3d4e5f6...`)
- **진짜 intermediate 파일들이 여기에 생성됨**

#### 구조
```
work/
├── a1/b2c3d4.../          # ALIGN_AND_SORT task
│   ├── sample.bam         # 실제 BAM 파일 생성
│   ├── sample.bam.bai
│   └── .command.sh        # 실행 스크립트
├── c4/d5e6f7.../          # MARK_DUPLICATES task
│   ├── sample.md.bam      # 실제 마크된 BAM
│   └── duplicate_metrics.txt
├── e8/f9g0h1.../          # CALL_VARIANTS task
│   └── sample_filtered.vcf.gz
└── ...
```

#### 중요한 점
- **모든 실제 계산과 파일 생성은 여기서 발생**
- publishDir는 이 work/ 파일들을 results/로 **복사 또는 링크**
- --cleanup 없이는 **무한정 쌓임** (디스크 풀 위험)

---

### 2단계: results/ 디렉토리 (publishDir)

#### publishDir의 두 가지 모드

##### A. mode: 'symlink' (심볼릭 링크)
```groovy
process ALIGN_AND_SORT {
    publishDir "${params.outdir}/alignment", mode: 'symlink'
    ...
}
```

**동작:**
```bash
results/alignment/sample.bam -> ../../work/a1/b2c3.../sample.bam
```

**특징:**
- ✅ 디스크 공간 절약 (링크만 생성)
- ⚠️ work/ 삭제 시 **파일 접근 불가** (broken symlink)
- 🎯 **대용량 파일에 사용** (BAM, 중간 결과)

**현재 사용:**
- `ALIGN_AND_SORT` (alignment/*.bam)
- `MARK_DUPLICATES` (alignment/*.md.bam)

##### B. mode: 'copy' (실제 복사)
```groovy
process CALL_VARIANTS {
    publishDir "${params.outdir}/variant", mode: 'copy'
    ...
}
```

**동작:**
```bash
cp work/e8/f9g0.../sample.vcf.gz results/variant/
```

**특징:**
- ✅ 독립적 파일 (work/ 삭제해도 안전)
- ❌ 디스크 공간 2배 사용
- 🎯 **최종 결과물에 사용** (VCF, 리포트)

**현재 사용:**
- `CALL_VARIANTS` (variant/)
- `MANTA_SV` (sv/)
- `GCNV_*` (cnv/)
- `PARAPHASE_RUN` (pseudogene/)
- `GENERATE_SUMMARY_REPORT` (summary/)
- 모든 분석 결과

---

### 3단계: final_archived_results/ (수동 아카이브)

#### 생성 방법
```bash
./scripts/archive_results.sh
```

#### 수행 작업
```bash
# BAM은 symlink를 실제 파일로 복사 (-L 옵션)
cp -L results/alignment/*.md.bam final_archived_results/bam/

# 나머지는 일반 복사
cp -r results/variant/* final_archived_results/vcf/
cp -r results/summary/* final_archived_results/summary/
...
```

#### 목적
- work/ 삭제 전에 **symlink를 실제 파일로 변환**
- 영구 보관용 최종 아카이브

---

## 🔄 전체 흐름 예시 (ALIGN_AND_SORT)

### 실행 과정
```bash
1. [work/] Nextflow가 task 실행
   work/a1/b2c3d4.../
   └── sample.bam (4GB 실제 파일 생성)

2. [results/] publishDir로 symlink 생성
   results/alignment/sample.bam -> ../../work/a1/b2c3.../sample.bam

3. [final_archived_results/] 아카이브 실행
   cp -L results/alignment/sample.bam final_archived_results/bam/
   → 4GB 실제 파일 복사됨
```

### 디스크 사용량
```
work/                    : 4GB (실제 파일)
results/alignment/       : 0GB (symlink만)
final_archived_results/  : 4GB (실제 파일 복사)
-------------------------------------------
총합                     : 8GB (중복)
```

### --cleanup 후
```
work/                    : 삭제됨
results/alignment/       : ⚠️ broken symlink (접근 불가)
final_archived_results/  : 4GB (유일한 실제 파일)
-------------------------------------------
총합                     : 4GB
```

---

## 🎯 문제점 및 해결 방안

### 현재 구조의 문제

#### 문제 1: work/ 디렉토리 크기 폭발
```
work/
├── a1/b2c3.../sample.bam (4GB)
├── c4/d5e6.../sample.md.bam (4GB)
├── e8/f9g0.../vcf 파일들 (500MB)
└── ... (수백 GB 누적)
```
→ 12 샘플 배치 시 **work/ 디렉토리만 180GB+**

#### 문제 2: 중복 저장
```
work/                    : 실제 파일 (180GB)
results/                 : symlink (0GB) + copy (20GB)
final_archived_results/  : 실제 파일 복사 (40GB)
-------------------------------------------
총합                     : 240GB
```

#### 문제 3: final_archived_results의 정체성 혼란
- 현재: "최종 아카이브"처럼 보임
- 실제: **results/를 단순 복사한 것** (일부는 symlink 해제)
- 혼란: intermediate vs final 구분 불명확

---

## ✅ 새로운 구조 제안

### 개선된 3단계 구조

```
[1단계] work/ → analysis/<work_dir>/<sample_name>/
   - Nextflow work dir는 그대로 사용
   - 중요 intermediate 파일만 선별 복사
   - 완료 후 work/ 삭제 가능

[2단계] analysis/ (모든 intermediate 파일)
   - publishDir의 기본 대상
   - symlink → copy로 변경 (독립성 확보)
   - work/ 삭제해도 안전

[3단계] output/ (포털용 최종 결과만)
   - summary/, snapshots/
   - 필터링된 VCF, 최종 BAM
   - 리포트, 시각화
```

### publishDir 모드 변경 제안

```groovy
// 기존 (문제)
publishDir "${params.outdir}/alignment", mode: 'symlink'

// 개선 (해결)
publishDir "${params.outdir}/alignment", mode: 'copy'
```

**장점:**
- ✅ work/ 삭제 가능 (디스크 확보)
- ✅ analysis/만 보존하면 재현 가능
- ✅ intermediate vs final 명확 분리

**단점:**
- ❌ 디스크 사용량 증가 (하지만 work/ 삭제로 상쇄)

---

## 📊 디스크 사용량 비교

### 현재 구조 (--cleanup 전)
```
work/               : 180GB (실제 작업)
results/            : 20GB (copy 파일만)
final_archived_results/ : 40GB (아카이브)
-------------------------------------------
총합                : 240GB
```

### 현재 구조 (--cleanup 후)
```
work/               : 삭제됨
results/            : ⚠️ symlink 깨짐
final_archived_results/ : 40GB (유일한 보존)
-------------------------------------------
총합                : 40GB (하지만 재분석 불가)
```

### 새로운 구조 (개선)
```
work/               : 180GB (작업 중) → 삭제
analysis/           : 60GB (중요 intermediate)
output/             : 10GB (포털용 최종)
-------------------------------------------
총합                : 70GB (재분석 가능)
```

---

## 🔧 구현 필요 사항

### 1. publishDir 모드 변경

#### align.nf
```groovy
// 변경 전
publishDir "${params.outdir}/alignment", mode: 'symlink'

// 변경 후
publishDir "${params.outdir}/alignment", mode: 'copy'
```

### 2. main.nf의 onComplete 수정

```groovy
workflow.onComplete {
    if (workflow.success) {
        // 1. output/ 복사 (이미 구현됨)
        copyToOutput()
        
        // 2. work/ 삭제 (--cleanup 시)
        if (params.cleanup) {
            deleteWorkDir()
        }
        
        // 3. final_archived_results 제거 (불필요)
        // archive_results.sh 사용 중단
    }
}
```

### 3. 불필요한 파일 제거

```bash
# 더 이상 필요 없음
rm scripts/archive_results.sh

# 기존 데이터 정리
rm -rf final_archived_results/  # analysis/로 통합
```

---

## 🎯 결론

### intermediate 파일의 진짜 위치

| 파일 종류 | 원본 위치 | 현재 저장 | 새 구조 |
|----------|----------|----------|---------|
| **실제 작업** | work/ | work/ (삭제됨) | work/ → 삭제 |
| **중간 결과** | work/ | results/ (symlink) | analysis/ (copy) |
| **최종 결과** | work/ | results/ (copy) | output/ (copy) |
| **아카이브** | - | final_archived_results/ | ❌ 불필요 |

### 핵심 깨달음

1. **final_archived_results는 intermediate가 아님**
   - results/를 단순 복사한 "백업"
   - work/ symlink를 실제 파일로 변환한 것
   
2. **진짜 intermediate는 work/에 있음**
   - 모든 계산과 파일 생성의 원본
   - publishDir는 이를 복사/링크만 함
   
3. **새 구조의 명확성**
   - `analysis/`: 모든 중간 파일 (재분석 가능)
   - `output/`: 포털용 최종 결과만
   - `work/`: 작업 완료 후 삭제 가능

