"""分布函数的时间演化与速度网格分辨率 (读 out/f_*.vpd；均匀与非均匀 v 网格均可)。

子命令:
  fv   RUN [REF] --d 0.25,20,500 --sp e,i [--t ...]
        给定离壁距离处 f(v) 随时间的演化：同一张图里各时刻叠画 (颜色 = 时间，蓝色由浅到深)，
        行 = 粒子，列 = 离壁距离；每张图底部的短竖线为该算例的速度格面 (看分布变窄时格是否够密)。
        给出 REF 时 REF 以虚线画在同一时刻 (速度网格可以不同)。
  res  RUN... [--sp e,i]
        分辨率监测：每个快照计算分布宽度 σ 与分布中心 u 处的速度格宽 ΔU 之比 σ/ΔU，
        画体区 (n > 1e-3 n_max 的格) 的最小值与贴壁格的值随时间的变化；虚线为 σ/Δv = 1.7 (电子 T 误差约 4%) 与 3.6 (约 1%)。
  phase RUN [--t 0.2,0.6,1,2,3,4] [--xrange 0,1]
        相空间 f(x, v)，x ∈ [0, L] (缺省中心到右壁)，画出每个相空间格的边 (非均匀 x、v 网格一目了然)；
        行 = 粒子，列 = 时刻，各图自己的色标 (对数，下限 --floor)，速度范围按该图分布自动取
  grid FILE_OR_RUN... [--ref runs/Case83]
        画速度网格的格宽 ΔU(U) (格面表文件或算例目录，热速度单位)；给出 --ref 时叠画参考算例的需求包络 σ_min(U)/k (k = 2、3)。
通用选项: --side avg|left|right (缺省 avg，两壁对折平均，v 取指向壁面为正)，--floor (对数纵轴下限，相对最大值，缺省 1e-6)，
          --vlim (速度范围；缺省自动)，--out，--show
速度以本粒子热速度为单位 (电子 v_Te，离子 v_Ti)。
例: python plot/VPevol.py fv runs/Case96 --d 0.25,20,500
    python plot/VPevol.py fv runs/Case96 runs/Case83 --d 0.25 --t 0.4,1,2,3,4
    python plot/VPevol.py res runs/Case95 runs/Case96 runs/Case97 runs/Case84
    python plot/VPevol.py grid case/vfaces/e128.txt case/vfaces/ei256.txt runs/Case84 --ref runs/Case83
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vpio import Run, moments, v_edges  # noqa: E402

import matplotlib  # noqa: E402

if "--show" not in sys.argv:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.cm import ScalarMappable  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, Normalize  # noqa: E402

from VPplot import CAT, DASH  # noqa: E402
from VPdist import cell_avg, mass_of, phase, rname  # noqa: E402

GRAY = "#8a8984"
TCMAP = LinearSegmentedColormap.from_list("vp_time", ["#b7d3f6", "#5598e7", "#2a78d6", "#184f95", "#0b2a55"])
SPNAME = {"e": "electrons", "i": "ions"}


def save(fig, args, default):
    out = args.out or default
    fig.savefig(out, dpi=args.dpi, bbox_inches="tight")
    print("图已保存:", out)
    if args.show:
        plt.show()


def snap_indices(run, tstr):
    if not tstr:
        return list(range(len(run.snap_files)))
    out = []
    for t in (float(v) for v in tstr.split(",")):
        k = int(np.argmin(np.abs(run.snap_t - t)))
        if abs(run.snap_t[k] - t) > 0.01:
            print(f"{rname(run)}: 没有 t={t:g} 的快照，用最近的 t={run.snap_t[k]:.4g}")
        out.append(k)
    return sorted(set(out))


def species_list(run, sps):
    s = run.snap(0)
    return [sp for sp in sps if f"f_{sp}" in s]


# ======================= fv: 给定位置 f(v) 的时间演化 =======================
def cmd_fv(args):
    run = Run(args.runs[0])
    ref = Run(args.runs[1]) if len(args.runs) > 1 else None
    sps = species_list(run, args.sp.split(","))
    ds = [float(v) for v in (args.d or "0.25,20,500").split(",")]
    ks = snap_indices(run, args.t)
    ts = run.snap_t[ks]
    norm = Normalize(ts.min(), ts.max() if ts.max() > ts.min() else ts.min() + 1)
    fig, axs = plt.subplots(len(sps), len(ds), figsize=(4.1 * len(ds) + 1.6, 3.4 * len(sps)),
                            squeeze=False, constrained_layout=True)
    skipped = []
    for r, sp in enumerate(sps):
        fmax = {d: 0.0 for d in ds}
        vl = {d: 0.0 for d in ds}
        for k in ks:
            s = run.snap(k)
            t = float(s["t"])
            X, v, ve, F, _, vlab = phase(run, s, sp, args.side)
            have_ref = False
            if ref is not None:
                s0 = ref.snap(t=t)
                have_ref = abs(float(s0["t"]) - t) <= 0.01  # 参考算例没有该时刻的快照时不画
                if have_ref:
                    X0, v0, ve0, F0, _, _ = phase(ref, s0, sp, args.side)
                elif r == 0:
                    skipped.append(t)
            col = TCMAP(norm(t))
            for c, d in enumerate(ds):
                ax = axs[r][c]
                i = int(np.argmin(np.abs(0.5 * (X[1:] + X[:-1]) - d)))
                fv = F[i]
                ax.plot(v, fv, color=col, lw=1.1)
                if have_ref:
                    fr = cell_avg(X0, F0, X[i], X[i + 1])
                    ax.plot(v0, fr, color=col, lw=0.9, ls=(0, (3, 2)))
                fmax[d] = max(fmax[d], fv.max())
                hit = np.abs(v[fv > args.floor * fv.max()]) if fv.max() > 0 else np.array([0.0])
                vl[d] = max(vl[d], hit.max())
                if k == ks[-1]:
                    ax.set_title(f"{SPNAME[sp]}, d = {0.5 * (X[i] + X[i + 1]):.3g} λ", fontsize=9)
        for c, d in enumerate(ds):
            ax = axs[r][c]
            lim = args.vlim or 1.15 * vl[d]
            ax.set_xlim(-lim, lim)
            if fmax[d] <= 0:  # 这一位置还没有等离子体
                ax.text(0.5, 0.5, "no plasma yet", transform=ax.transAxes, ha="center", color=GRAY)
                continue
            if not args.lin:
                ax.set_yscale("log")
                ax.set_ylim(args.floor * fmax[d], 2.5 * fmax[d])
                y0 = args.floor * fmax[d]
                ytick = (y0, y0 * 3)
            else:
                ytick = (0, 0.04 * fmax[d])
            # 速度格面 (本算例)：底部短竖线
            ve = phase(run, run.snap(ks[-1]), sp, args.side)[2]
            sel = ve[(ve >= -lim) & (ve <= lim)]
            ax.vlines(sel, *ytick, color=GRAY, lw=0.4)
            ax.set_xlabel(vlab)
            if c == 0:
                ax.set_ylabel(f"f_{sp}")
    if skipped:
        print(f"{rname(ref)} 没有这些时刻的快照，未画参考: " + ", ".join(f"{t:.3g}" for t in skipped))
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=TCMAP), ax=axs, shrink=0.8, location="right")
    cb.set_label("t [τi]")
    title = f"{rname(run)}: f(v) evolution (side={args.side}); ticks at bottom = velocity cell faces"
    if ref is not None:
        title += f"; dashed = {rname(ref)} (avg over the same x cell)"
    fig.suptitle(title, fontsize=10)
    save(fig, args, os.path.join(run.path, "evol_fv.png"))


# ======================= res: σ/ΔU 随时间 =======================
def resolution(run, sps):
    """返回 {sp: (t, 体区最小 σ/ΔU, 贴壁格 σ/ΔU, 体区最小 σ)}，速度以本粒子热速度为单位。"""
    out = {sp: ([], [], [], []) for sp in sps}
    for k in range(len(run.snap_files)):
        s = run.snap(k)
        for sp in sps:
            if f"f_{sp}" not in s:
                continue
            m = mass_of(s, sp)
            sc = np.sqrt(m) if sp == "i" else 1.0
            f = np.asarray(s[f"f_{sp}"], float) / sc
            v = np.asarray(s[f"v_{sp}"]) * sc
            ve = v_edges(s, sp) * sc
            dv = np.diff(ve)
            n, u, T, _ = moments(f, v, 1.0, dv=dv)
            sig = np.sqrt(np.maximum(T, 0.0))
            j = np.clip(np.searchsorted(ve, u) - 1, 0, len(dv) - 1)
            ratio = sig / dv[j]
            ok = n > 1e-3 * n.max()
            if not ok.any():
                continue
            o = out[sp]
            o[0].append(float(s["t"]))
            o[1].append(ratio[ok].min())
            o[2].append(0.5 * (ratio[0] + ratio[-1]))
            o[3].append(sig[ok].min())
    return {sp: tuple(np.array(a) for a in o) for sp, o in out.items()}


def cmd_res(args):
    runs = [Run(p) for p in args.runs]
    sps = args.sp.split(",")
    fig, axs = plt.subplots(1, len(sps), figsize=(5.2 * len(sps), 3.8), squeeze=False, constrained_layout=True)
    print(f"{'算例':10s} {'粒子':4s} " + " ".join(f"{'t=' + str(t):>12s}" for t in (0.5, 1, 2, 3, 4)) + "   (体区最小 σ/ΔU)")
    for ri, run in enumerate(runs):
        res = resolution(run, sps)
        for c, sp in enumerate(sps):
            t, rb, rw, _ = res[sp]
            if not len(t):
                continue
            ax = axs[0][c]
            ax.plot(t, rb, color=CAT[ri % len(CAT)], lw=1.6, marker="o", ms=3, label=f"{rname(run)} bulk min")
            ax.plot(t, rw, color=CAT[ri % len(CAT)], lw=1.0, ls=(0, (3, 2)), label=f"{rname(run)} wall cell")
            row = []
            for tt in (0.5, 1, 2, 3, 4):
                kk = int(np.argmin(np.abs(t - tt)))
                row.append(f"{rb[kk]:12.2f}" if abs(t[kk] - tt) < 0.15 else f"{'':>12}")
            print(f"{rname(run):10s} {sp:4s} " + " ".join(row))
    for c, sp in enumerate(sps):
        ax = axs[0][c]
        for y, lab in ((1.7, "σ/Δv = 1.7"), (3.6, "3.6")):
            ax.axhline(y, color=GRAY, lw=0.8, ls=":")
            ax.text(0.01, y, lab, transform=ax.get_yaxis_transform(), fontsize=7, color=GRAY, va="bottom")
        ax.set_yscale("log")
        ax.set_xlabel("t [τi]")
        ax.set_ylabel("σ / ΔU at distribution centre")
        ax.set_title(SPNAME[sp])
        ax.legend(fontsize=7)
    save(fig, args, os.path.join(runs[0].path, "evol_res.png"))


# ======================= grid: ΔU(U) =======================
def load_faces(item):
    """格面表文件或算例目录 → (名字, {sp: 格面 (热速度单位)})。"""
    if os.path.isfile(item):
        f = []
        for line in open(item):
            line = line.split("#")[0]
            f += [float(a) for a in line.split()]
        return os.path.basename(item), {"e": np.array(f)}
    run = Run(item)
    s = run.snap(0)
    out = {}
    for sp in ("e", "i"):
        if f"v_{sp}" in s:
            sc = np.sqrt(mass_of(s, sp)) if sp == "i" else 1.0
            out[sp] = v_edges(s, sp) * sc
    return rname(run), out


def cmd_grid(args):
    fig, ax = plt.subplots(figsize=(6.6, 4.2), constrained_layout=True)
    ci = 0
    for item in args.runs:
        name, faces = load_faces(item)
        for sp, fc in faces.items():
            w = np.diff(fc)
            c = 0.5 * (fc[1:] + fc[:-1])
            m = c >= 0
            lab = f"{name}" + (f" ({sp})" if len(faces) > 1 else "") + f", N={len(w)}"
            ax.step(c[m], w[m], where="mid", color=CAT[ci % len(CAT)], lw=1.4,
                    ls=DASH[0] if sp == "e" else DASH[1], label=lab)
            ci += 1
    if args.ref:
        sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tools"))
        from make_vfaces import envelope
        for sp, col in (("e", "#9a9893"), ("i", "#3d3c39")):  # 需求包络用灰色系，与网格曲线区分
            U, req = envelope(args.ref, [sp], 1.0, 0.15, 6.0, ngrid=3001)
            for kf, ls in ((2, "-"), (3, ":")):
                ax.plot(U, req / kf, color=col, lw=1.0, ls=ls,
                        label=f"need σ/{kf} ({SPNAME[sp]}, {os.path.basename(args.ref.rstrip('/'))})")
    ax.set_yscale("log")
    ax.set_xlabel("|U| [thermal velocity of each species]")
    ax.set_ylabel("ΔU (cell width)")
    ax.set_xlim(0, 6)
    ax.legend(fontsize=7)
    save(fig, args, args.out or "vgrid.png")


# ======================= phase: 相空间 f(x, v) 与网格 =======================
def cmd_phase(args):
    """x ∈ [x0, x1] (L 单位) 的相空间 f(x, v)，画出每个相空间格的边 (看网格)；行 = 粒子，列 = 时刻。"""
    from matplotlib.colors import LogNorm

    from VPplot import SEQ
    run = Run(args.runs[0])
    sps = species_list(run, args.sp.split(","))
    ks = snap_indices(run, args.t or "0.2,0.6,1,2,3,4")
    x0, x1 = (float(a) for a in args.xrange.split(","))
    fig, axs = plt.subplots(len(sps), len(ks), figsize=(3.3 * len(ks) + 0.8, 3.0 * len(sps)),
                            squeeze=False, constrained_layout=True)
    for r, sp in enumerate(sps):
        for c, k in enumerate(ks):
            s = run.snap(k)
            m = mass_of(s, sp)
            sc = np.sqrt(m) if sp == "i" else 1.0
            f = np.asarray(s[f"f_{sp}"], float) / sc
            ve = v_edges(s, sp) * sc
            xf = np.asarray(s["x_face"])
            i0 = max(int(np.searchsorted(xf, x0 + 1e-12)) - 1, 0)
            i1 = min(int(np.searchsorted(xf, x1 - 1e-12)), len(xf) - 1)
            X, F = xf[i0:i1 + 1], f[i0:i1]
            fmax = F.max()
            norm = LogNorm(vmin=args.floor * fmax, vmax=fmax)
            ax = axs[r][c]
            ax.pcolormesh(X, ve, np.maximum(F.T, args.floor * fmax), norm=norm, cmap=SEQ, shading="flat",
                          edgecolors=(0.3, 0.3, 0.3, 0.35) if not args.nogrid else "none",
                          linewidth=0.15, rasterized=True)
            hit = np.abs(0.5 * (ve[1:] + ve[:-1])[(F > args.floor * fmax).any(axis=0)])
            vl = args.vlim or (1.15 * hit.max() if hit.size else np.abs(ve).max())
            ax.set_ylim(-vl, vl)
            ax.set_xlim(X[0], X[-1])
            ax.grid(False)
            ax.set_title(f"{SPNAME[sp]}, t = {float(s['t']):.3g} τi", fontsize=9)
            if r == len(sps) - 1:
                ax.set_xlabel("x [L]")
            if c == 0:
                ax.set_ylabel("v [v_Te]" if sp == "e" else "v [v_Ti]")
            if c == len(ks) - 1:
                fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=SEQ), ax=ax).set_label(f"f_{sp} (each panel own scale)")
    fig.suptitle(f"{rname(run)}: f(x, v) for x ∈ [{x0:g}, {x1:g}] L, cell edges drawn "
                 f"(colour floor {args.floor:g} of panel max)", fontsize=10)
    save(fig, args, os.path.join(run.path, "evol_phase.png"))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["fv", "res", "grid", "phase"])
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--t", help="fv: 时刻 (τi)，逗号分隔；缺省为全部快照")
    ap.add_argument("--sp", default="e,i")
    ap.add_argument("--d", help="fv: 离壁距离 (λ)，逗号分隔，缺省 0.25,20,500")
    ap.add_argument("--side", default="avg", choices=["avg", "left", "right"])
    ap.add_argument("--ref", help="grid: 参考算例 (画需求包络)")
    ap.add_argument("--xrange", default="0,1", help="phase: x 范围 (L)，缺省 0,1 (中心到右壁)")
    ap.add_argument("--nogrid", action="store_true", help="phase: 不画格边")
    ap.add_argument("--lin", action="store_true", help="fv: 线性纵轴")
    ap.add_argument("--floor", type=float, default=1e-6)
    ap.add_argument("--vlim", type=float)
    ap.add_argument("--out")
    ap.add_argument("--dpi", type=int, default=150)
    ap.add_argument("--show", action="store_true")
    args = ap.parse_args()
    {"fv": cmd_fv, "res": cmd_res, "grid": cmd_grid, "phase": cmd_phase}[args.cmd](args)


if __name__ == "__main__":
    main()
