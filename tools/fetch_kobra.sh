#!/bin/bash
# 把 kobra 的诊断输出拷贝到 VPcode/runs/kobra_<标签>/ 下（只拷贝，不删除 kobra 里的原文件）
# 用法: tools/fetch_kobra.sh <标签> [kobra 数据目录] [kobra 算例 .dat 文件]
#   数据目录缺省为 $kobra_data，否则 ~/Vlasov/quasi-kobra-sup/data
# 例:   tools/fetch_kobra.sh n1024 ~/Vlasov/quasi-kobra-sup/data ~/Vlasov/quasi-kobra-sup/naga/dts.dat
set -e
TAG=$1
if [ -z "$TAG" ]; then
    echo "用法: tools/fetch_kobra.sh <标签> [kobra 数据目录] [kobra 算例 .dat 文件]"
    exit 1
fi
SRC=${2:-${kobra_data:-$HOME/Vlasov/quasi-kobra-sup/data}}
ROOT=$(cd "$(dirname "$0")/.." && pwd)
DST=$ROOT/runs/kobra_$TAG
if [ ! -f "$SRC/diag_budget.dat" ]; then
    echo "在 $SRC 下没有 diag_budget.dat"
    exit 1
fi
mkdir -p "$DST"
for f in "$SRC"/diag_*.dat "$SRC"/dtsmon.dat "$SRC"/phimon_sum.dat; do
    [ -f "$f" ] && cp -p "$f" "$DST/" && echo "拷贝 $(basename "$f")"
done
if [ -n "$3" ]; then
    cp -p "$3" "$DST/" && echo "拷贝算例文件 $(basename "$3")"
fi
echo "完成: $DST"
echo "画图: ~/plotenv/bin/python plot/KobraDiag.py $DST [runs/Case25]"
