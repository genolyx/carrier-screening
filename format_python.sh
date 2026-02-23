#!/bin/bash
# Python 코드 자동 포맷팅 스크립트

FILE="${1:-dashboard/app.py}"

echo "Formatting $FILE with autopep8..."

# autopep8 설치 확인
if ! command -v autopep8 &> /dev/null; then
    echo "Installing autopep8..."
    pip3 install --user autopep8
fi

# 백업 생성
cp "$FILE" "${FILE}.backup_$(date +%Y%m%d_%H%M%S)"

# 자동 포맷팅
autopep8 --in-place --aggressive --aggressive "$FILE"

echo "✅ Formatted! Backup created."
