"""壁面附近的网格与分布函数 (AMR / 非均匀 x 网格算例；读 RUN/out/f_<ts>.vpd)。

用法:
  python plot/KobraWallGrid.py RUN [--t 0.2,0.8,1.4] [--dmax 64] [--side left|right] [--out 前缀] [--show]
例:
  python plot/KobraWallGrid.py ~/Vlasov/kobra-amr/naga --t 0.2,0.8,1.4

输出两张图 (缺省存在 RUN 下，前缀 wallgrid):
  _phase.png  第 1 行: 格宽 Δx 随离壁距离的分布 (按格宽着色)；第 2、3 行: 各时刻电子、离子的相空间 f(d, v)，
              按实际 (非均匀) 格画出，竖线为格面，灰色虚线为粗细界面
  _fv.png     最后一个时刻，每个粗细界面两侧各一个格的 f(v) (电子、离子，对数轴)，看分布在界面处是否连续
坐标: d 为离壁距离 (λ)；v 取指向壁面为正；电子 v 以 v_Te 计，离子以 v_Ti = v_Te/√(m_i/m_e) 计 (f 相应换算，∫f dv = n)。
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vpio import VPD, Run  # noqa: E402

import matplotlib  # noqa: E402

if "--show" not in sys.argv:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

from VPplot import CAT, SEQ  # noqa: E402

MI = 1836.0


def side_view(s, side, lam):
    """返回一侧壁面的格面 (离壁距离, λ)、格宽、f_e/f_i (行 = 离壁由近到远，列 = 指向壁面为正的 v)。"""
    xf = np.asarray(s["x_face"])
    ve = np.asarray(s["v_e"])
    vi = np.asarray(s["v_i"]) * np.sqrt(MI)
    nx, nv = len(xf) - 1, len(ve)
    fe = np.asarray(s["f_e"]).reshape(nx, nv)
    fi = np.asarray(s["f_i"]).reshape(nx, nv) / np.sqrt(MI)
    if side == "left":
        df = (xf + 1.0) / lam
        # 左壁：指向壁面为 −v，速度轴翻转
        fe, fi = fe[:, ::-1], fi[:, ::-1]
        ve, vi = -ve[::-1], -vi[::-1]
    else:
        df = ((1.0 - xf) / lam)[::-1]
        fe, fi = fe[::-1], fi[::-1]
    return df, np.diff(df), ve, vi, fe, fi


def vfaces(v):
    dv = v[1] - v[0]
    return np.concatenate([v - 0.5 * dv, [v[-1] + 0.5 * dv]])


def vrange(v, f, floor=1e-6):
    m = (f > floor * f.max()).any(axis=0)
    if not m.any():
        return v[0], v[-1]
    lo, hi = v[m][0], v[m][-1]
    pad = 0.1 * (hi - lo)
    return lo - pad, hi + pad


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run")
    ap.add_argument("--t", default=None, help="时刻 (τi)，逗号分隔；缺省取 3 个快照 (第一个非零、中间、最后)")
    ap.add_argument("--dmax", type=float, default=64.0, help="离壁距离上限 (λ)")
    ap.add_argument("--side", default="left", choices=["left", "right"])
    ap.add_argument("--floor", type=float, default=1e-6, help="对数色标下限 (相对最大值)")
    ap.add_argument("--out", default=None)
    ap.add_argument("--show", action="store_true")
    a = ap.parse_args()

    run = Run(os.path.expanduser(a.run))
    ts = np.asarray(run.snap_t)
    if a.t:
        want = [float(x) for x in a.t.split(",")]
    else:
        nz = ts[ts > 0]
        want = [nz[0], nz[len(nz) // 2], nz[-1]] if len(nz) >= 3 else list(nz)
    snaps = [run.snap(t=t) for t in want]
    lam = float(VPD(run.snap_files[0]).meta.get("lambda", 1e-3))

    # ---- 图 1: 格宽 + 相空间 ----
    nt = len(snaps)
    fig, ax = plt.subplots(3, nt, figsize=(5.2 * nt, 11), gridspec_kw={"height_ratios": [0.7, 1, 1]})
    ax = np.atleast_2d(ax).reshape(3, nt)
    for j, s in enumerate(snaps):
        df, dx, ve, vi, fe, fi = side_view(s, a.side, lam)
        m = df[1:] <= a.dmax + 1e-9
        k = np.where(m)[0]
        dfm = df[: k[-1] + 2]
        widths = np.unique(np.round(dx[k], 6))
        cmap = {w: CAT[i % len(CAT)] for i, w in enumerate(widths)}
        jumps = [df[i + 1] for i in k[:-1] if abs(dx[i + 1] - dx[i]) > 1e-6]
        # 第 1 行：格宽
        g = ax[0, j]
        for i in k:
            g.bar(0.5 * (df[i] + df[i + 1]), dx[i], width=dx[i], color=cmap[round(dx[i], 6)],
                  edgecolor="white", linewidth=0.6)
        handles = [Patch(color=cmap[w], label=f"Δx = {w:g} λ") for w in widths]
        for xj in jumps:
            g.axvline(xj, color="#8a8984", ls="--", lw=0.8)
        g.set_xlim(0, a.dmax)
        g.set_ylabel("Δx [λ]")
        g.set_title(f"t = {float(np.asarray(s['t'])):.2f} τi: cell widths ({len(k)} cells)", fontsize=10)
        if j == 0:
            g.legend(handles=handles, frameon=False, fontsize=8)
        # 第 2、3 行：相空间
        for r, (v, f, sp, vu) in enumerate([(ve, fe, "f_e", "v_Te"), (vi, fi, "f_i", "v_Ti")], start=1):
            g = ax[r, j]
            fm = f[: len(dfm) - 1]
            vmax = fm.max()
            pc = g.pcolormesh(dfm, vfaces(v), np.clip(fm.T, a.floor * vmax, None), cmap=SEQ,
                              norm=LogNorm(vmin=a.floor * vmax, vmax=vmax), shading="flat")
            for xf in dfm:
                g.axvline(xf, color="white", lw=0.35, alpha=0.8)
            for xj in jumps:
                g.axvline(xj, color="#8a8984", ls="--", lw=0.8)
            g.set_ylim(*vrange(v, fm, a.floor))
            g.axhline(0, color="#8a8984", lw=0.5)
            g.set_xlim(0, a.dmax)
            g.set_ylabel(f"v (toward wall +) [{vu}]")
            g.set_title(sp, fontsize=10)
            fig.colorbar(pc, ax=g, pad=0.01)
            if r == 2:
                g.set_xlabel("distance from wall [λ]")
    fig.suptitle(f"{os.path.basename(os.path.normpath(a.run))}: grid and phase-space f near the {a.side} wall (white: cell faces, dashed: coarse/fine interfaces)")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    pre = a.out or os.path.join(os.path.expanduser(a.run), "wallgrid")
    fig.savefig(pre + "_phase.png", dpi=130)
    print(f"图已保存: {pre}_phase.png")

    # ---- 图 2: 界面两侧的 f(v) ----
    s = snaps[-1]
    df, dx, ve, vi, fe, fi = side_view(s, a.side, lam)
    k = np.where(df[1:] <= a.dmax + 1e-9)[0]
    pairs = [(i, i + 1) for i in k[:-1] if abs(dx[i + 1] - dx[i]) > 1e-6]
    fig, ax = plt.subplots(2, max(len(pairs), 1), figsize=(4.2 * max(len(pairs), 1), 7), squeeze=False)
    for j, (i0, i1) in enumerate(pairs):
        for r, (v, f, sp, vu) in enumerate([(ve, fe, "f_e", "v_Te"), (vi, fi, "f_i", "v_Ti")]):
            g = ax[r, j]
            for c, i in zip([CAT[0], CAT[1]], [i0, i1]):
                dc = 0.5 * (df[i] + df[i + 1])
                g.semilogy(v, np.maximum(f[i], 1e-30), color=c, lw=1.3, label=f"d = {dc:g} λ (Δx {dx[i]:g})")
            lo, hi = vrange(v, f[[i0, i1]], a.floor)
            g.set_xlim(lo, hi)
            top = max(f[i0].max(), f[i1].max())
            g.set_ylim(a.floor * top, 2 * top)
            # 两侧密度与温度的相对差
            dv = v[1] - v[0]
            n0, n1 = f[i0].sum() * dv, f[i1].sum() * dv
            u0, u1 = (f[i0] * v).sum() * dv / n0, (f[i1] * v).sum() * dv / n1
            T0 = (f[i0] * (v - u0) ** 2).sum() * dv / n0
            T1 = (f[i1] * (v - u1) ** 2).sum() * dv / n1
            g.set_title(f"{sp}: interface d = {df[i1]:g} λ\nΔn/n = {(n1 / n0 - 1) * 100:+.2f}%, ΔT/T = {(T1 / T0 - 1) * 100:+.2f}%",
                        fontsize=9)
            g.set_xlabel(f"v (toward wall +) [{vu}]")
            g.grid(alpha=0.3)
            g.legend(frameon=False, fontsize=7)
    fig.suptitle(f"t = {float(np.asarray(s['t'])):.2f} τi: f(v) in the two cells adjacent to each coarse/fine interface")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(pre + "_fv.png", dpi=130)
    print(f"图已保存: {pre}_fv.png")
    if a.show:
        plt.show()


if __name__ == "__main__":
    main()
