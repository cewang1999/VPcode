#!/bin/bash
# 把 kobra-amr 算例的 VPD 输出拷贝到 VPcode/runs/<NAME>/ 下（只拷贝，不删除 kobra-amr 里的原文件）
# 用法: tools/fetch_kobra_amr.sh <NAME> [源目录]
#   NAME    算例名 (如 KB16N256)，即 runs/ 下的目录名；不能以 Case 开头 (避免与 VPcode 算例混淆)
#   源目录  缺省按顺序找:
#             $KOBRA_AMR/runs/<NAME>   (amr_cases/run.sh 的运行目录)
#             $KOBRA_AMR/naga          (直接在 naga 下运行；out/ 中是最后一次运行的结果)
#           KOBRA_AMR 缺省为 ~/Vlasov/kobra-amr
# 拷贝: out/ (f_*.vpd、prof.vpd)、flux_mpi_nrel100.out、<NAME>.dat、log、time.txt、cmd (有则拷)
# 算例还在运行时也可以拷，之后再运行一次即更新 (未写完的 VPD 记录 / 快照会被 vpio 跳过)
# 例:   tools/fetch_kobra_amr.sh KB16N256
#       tools/fetch_kobra_amr.sh KA4 ~/Vlasov/kobra-amr/runs/KA4
set -e
NAME=$1
if [ -z "$NAME" ]; then
    echo "用法: tools/fetch_kobra_amr.sh <NAME> [源目录]"
    exit 1
fi
case "$NAME" in
    Case*) echo "NAME 不能以 Case 开头 (runs/Case* 是 VPcode 的算例)"; exit 1 ;;
esac
KOBRA=${KOBRA_AMR:-$HOME/Vlasov/kobra-amr}
if [ -n "$2" ]; then
    SRC=$2
elif [ -d "$KOBRA/runs/$NAME/out" ]; then
    SRC=$KOBRA/runs/$NAME
else
    SRC=$KOBRA/naga
fi
SRC=$(cd "$SRC" && pwd)
if [ ! -f "$SRC/out/prof.vpd" ]; then
    echo "在 $SRC/out 下没有 prof.vpd (算例需在 .dat 中加 vpdprof)"
    exit 1
fi

# 算例文件：源目录下的 <NAME>.dat，没有时用 amr_cases/<NAME>.dat
DAT=$SRC/$NAME.dat
[ -f "$DAT" ] || DAT=$KOBRA/amr_cases/$NAME.dat
[ -f "$DAT" ] || DAT=

meta() { head -c 4000 "$SRC/out/prof.vpd" | LC_ALL=C grep -a "^meta $1 " | head -1 | awk '{print $3}'; }

ROOT=$(cd "$(dirname "$0")/.." && pwd)
DST=$ROOT/runs/$NAME
mkdir -p "$DST/out"
echo "源: $SRC"
echo "    prof.vpd: created $(meta created)，nx $(meta nx)，nv $(meta nv)，levels $(meta levels)，lambda $(meta lambda)"
n=0
for f in "$SRC"/out/*.vpd; do
    cp -p "$f" "$DST/out/" && n=$((n + 1))
done
echo "拷贝 out/ 下 $n 个 VPD 文件"
for f in flux_mpi_nrel100.out log time.txt cmd; do
    [ -f "$SRC/$f" ] && cp -p "$SRC/$f" "$DST/" && echo "拷贝 $f"
done
if [ -n "$DAT" ]; then
    cp -p "$DAT" "$DST/" && echo "拷贝算例文件 $(basename "$DAT")"
fi
echo "完成: $DST ($(du -sh "$DST" | cut -f1))"
echo "画图: MPLBACKEND=Agg ~/plotenv/bin/python plot/KobraFig7.py runs/$NAME runs/Case83 runs/Case87"
echo "      ~/plotenv/bin/python plot/KobraFlux.py runs/$NAME"
echo "      ~/plotenv/bin/python plot/VPcmp.py runs/$NAME runs/Case83   (需双精度算例，格面才能对齐)"
