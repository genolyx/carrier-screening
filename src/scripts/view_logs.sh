#!/bin/bash
# Nextflow 로그 뷰어 (서버에서 직접 사용)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(dirname "$SCRIPT_DIR")"
LOG_BASE="$BASE_DIR/log"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

function usage() {
    cat << EOF
Usage: $(basename $0) [OPTIONS]

Nextflow 로그 조회 도구

OPTIONS:
    -l, --list              모든 로그 파일 목록 표시
    -w, --work-dir DIR      Work directory (예: 2601)
    -s, --sample NAME       Sample name
    -f, --follow            실시간 tail -f 모드
    -n, --lines N           마지막 N줄만 표시 (기본: 전체)
    -h, --help              이 도움말 표시

EXAMPLES:
    # 로그 목록 보기
    $(basename $0) --list
    
    # 특정 샘플 로그 전체 보기
    $(basename $0) -w 2601 -s Sample_A10
    
    # 마지막 100줄만 보기
    $(basename $0) -w 2601 -s Sample_A10 -n 100
    
    # 실시간 모니터링 (tail -f)
    $(basename $0) -w 2601 -s Sample_A10 -f
    
    # 최신 샘플 자동 선택하여 tail
    $(basename $0) --latest -f

EOF
    exit 0
}

function list_logs() {
    echo -e "${BLUE}=== Available Nextflow Logs ===${NC}\n"
    
    if [ ! -d "$LOG_BASE" ]; then
        echo -e "${RED}Log directory not found: $LOG_BASE${NC}"
        exit 1
    fi
    
    local count=0
    for work_dir in $(ls -1 "$LOG_BASE" 2>/dev/null | sort); do
        local work_path="$LOG_BASE/$work_dir"
        [ ! -d "$work_path" ] && continue
        
        for sample_name in $(ls -1 "$work_path" 2>/dev/null | sort); do
            local log_file="$work_path/$sample_name/nextflow.log"
            if [ -f "$log_file" ]; then
                count=$((count + 1))
                local size=$(du -h "$log_file" | cut -f1)
                local modified=$(stat -c '%y' "$log_file" 2>/dev/null || stat -f '%Sm' "$log_file" 2>/dev/null)
                local modified_short=$(echo "$modified" | cut -d' ' -f1,2 | cut -d'.' -f1)
                
                # 완료 여부 확인
                local status="${YELLOW}RUNNING${NC}"
                if grep -q "Execution complete" "$log_file" 2>/dev/null; then
                    if grep -q "Goodbye" "$log_file" 2>/dev/null; then
                        status="${GREEN}SUCCESS${NC}"
                    else
                        status="${RED}FAILED${NC}"
                    fi
                fi
                
                echo -e "${count}. ${BLUE}$work_dir/$sample_name${NC}"
                echo -e "   Status: $status | Size: $size | Modified: $modified_short"
                echo ""
            fi
        done
    done
    
    if [ $count -eq 0 ]; then
        echo -e "${YELLOW}No log files found${NC}"
    else
        echo -e "${GREEN}Total: $count log files${NC}"
    fi
}

function get_latest_log() {
    # 가장 최근에 수정된 로그 파일 찾기
    local latest=$(find "$LOG_BASE" -name "nextflow.log" -type f -printf '%T@ %p\n' 2>/dev/null | sort -rn | head -1 | cut -d' ' -f2)
    
    if [ -z "$latest" ]; then
        echo -e "${RED}No log files found${NC}"
        exit 1
    fi
    
    echo "$latest"
}

function show_log() {
    local work_dir=$1
    local sample_name=$2
    local follow=$3
    local lines=$4
    
    local log_file="$LOG_BASE/$work_dir/$sample_name/nextflow.log"
    
    if [ ! -f "$log_file" ]; then
        echo -e "${RED}Log file not found: $log_file${NC}"
        exit 1
    fi
    
    echo -e "${BLUE}=== Nextflow Log: $work_dir/$sample_name ===${NC}"
    echo -e "${YELLOW}File: $log_file${NC}"
    echo -e "${YELLOW}Size: $(du -h "$log_file" | cut -f1)${NC}"
    echo ""
    
    if [ "$follow" = "true" ]; then
        # 실시간 tail
        echo -e "${GREEN}[Following log in real-time... Press Ctrl+C to exit]${NC}\n"
        tail -f -n ${lines:-100} "$log_file"
    else
        # 일반 출력
        if [ -n "$lines" ]; then
            tail -n $lines "$log_file"
        else
            cat "$log_file"
        fi
    fi
}

# Parse arguments
LIST=false
WORK_DIR=""
SAMPLE_NAME=""
FOLLOW=false
LINES=""
LATEST=false

while [[ $# -gt 0 ]]; do
    case $1 in
        -l|--list)
            LIST=true
            shift
            ;;
        -w|--work-dir)
            WORK_DIR="$2"
            shift 2
            ;;
        -s|--sample)
            SAMPLE_NAME="$2"
            shift 2
            ;;
        -f|--follow)
            FOLLOW=true
            shift
            ;;
        -n|--lines)
            LINES="$2"
            shift 2
            ;;
        --latest)
            LATEST=true
            shift
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo -e "${RED}Unknown option: $1${NC}"
            usage
            ;;
    esac
done

# Execute
if [ "$LIST" = true ]; then
    list_logs
elif [ "$LATEST" = true ]; then
    LATEST_LOG=$(get_latest_log)
    # Extract work_dir and sample_name from path
    WORK_DIR=$(echo "$LATEST_LOG" | awk -F'/' '{print $(NF-2)}')
    SAMPLE_NAME=$(echo "$LATEST_LOG" | awk -F'/' '{print $(NF-1)}')
    echo -e "${GREEN}Latest log: $WORK_DIR/$SAMPLE_NAME${NC}\n"
    show_log "$WORK_DIR" "$SAMPLE_NAME" "$FOLLOW" "$LINES"
elif [ -n "$WORK_DIR" ] && [ -n "$SAMPLE_NAME" ]; then
    show_log "$WORK_DIR" "$SAMPLE_NAME" "$FOLLOW" "$LINES"
else
    echo -e "${RED}Error: Missing required arguments${NC}\n"
    usage
fi
