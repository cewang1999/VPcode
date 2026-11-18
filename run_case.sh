#!/bin/bash
# 用法: ./run_case.sh <case编号> [线程数]
# 在 runs/Case<编号>/ 下后台运行，输出 runtime_diags.csv、out/ (VPD)、log.txt
set -e
N=$1
T=${2:-8}
if [ -z "$N" ]; then
    echo "用法: ./run_case.sh <case编号> [线程数]"
    exit 1
fi
ROOT=$(cd "$(dirname "$0")" && pwd)
D=$ROOT/runs/Case$N
mkdir -p "$D/case"
cp "$ROOT/case/Case$N.dat" "$D/case/"
# 非均匀 v 网格的格面表 (VFacesE / VFacesI 以运行目录为基准的相对路径 case/vfaces/...)
if grep -qE '^VFaces[EI]' "$ROOT/case/Case$N.dat"; then
    cp -r "$ROOT/case/vfaces" "$D/case/"
fi
cp -r "$ROOT/input" "$D/"
cd "$D"
OMP_NUM_THREADS=$T nohup "$ROOT/vp" "$N" > log.txt 2>&1 &
echo "Case$N 已在后台运行 (PID $!, $T 线程): $D"
