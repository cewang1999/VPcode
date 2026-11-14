"""VPcode 分布函数可视化 (读 out/f_<ts>.vpd；旧 txt 算例只有电子)。

子命令 (RUN 为算例目录):
  panel RUN    --t 0.3,0.5,1 --sp e,i        相空间 f(x, v) 小多图: 行 = 粒子, 列 = 时刻 (同一粒子共用色标)
  anim  RUN    --sp e,i [--fps 4]            全部快照的相空间动画 (GIF)
  slice RUN... --t 0.5 --d 0.25,10,100       给定离壁距离处的 f(v)；--maxw 叠加同 n、u、T 的麦克斯韦分布
  vt    RUN    --d 0.25 --sp e               给定离壁距离处 f(v) 随时间的演化 (速度 × 时间图)
  diff  RUN RUN0 --t 0.5 --sp e              两个算例同一时刻 f 之差 (同一网格)
  wall  TEST REF --d 0.25 --sp e [--t ...]   离壁 d 处 f(v) 逐时刻对比 (上: f，下: 差)；REF 按 TEST 的格做格宽加权平均，可比不同网格
通用选项:
  --side avg|left|right|full  缺省 avg (两壁对折平均)；对折时 v 取指向壁面为正
  --dmax 离壁距离上限 (λ)，--logd 距离对数轴，--lin 线性色标 (缺省对数)，--floor 对数色标下限 (相对最大值, 缺省 1e-6)
  --vlim 速度范围 (电子 v_Te, 离子 v_Ti；缺省自动)，--out 文件名，--show 弹窗
离子速度以 v_Ti = v_Te/√(m_i/m_e) 为单位 (f 相应乘 √(m_i/m_e)⁻¹，∫f dv 仍为 n)。
例: python plot/VPdist.py panel runs/Case79 --t 0.2,0.5,1.0,1.6 --dmax 200 --logd
    python plot/VPdist.py slice runs/Case79 runs/Case63 --t 1.0 --d 0.5,20,200 --sp i --maxw
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vpio import Run, moments  # noqa: E402

import matplotlib  # noqa: E402

if "--show" not in sys.argv:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.animation import FuncAnimation, PillowWriter  # noqa: E402
from matplotlib.colors import LogNorm, Normalize, TwoSlopeNorm  # noqa: E402

from VPplot import BLUE, CAT, DASH, DIV, SEQ  # noqa: E402  (同一套配色与样式)

GRAY = "#8a8984"


def rname(run):
    return os.path.basename(run.path.rstrip("/"))


def lam_of(run, s):
    return s.meta.get("lambda", run.case.get("lambda", 1e-3))


def mass_of(s, sp):
    return float(s.meta.get(f"mass_{sp}", 1.0 if sp == "e" else 1836.0))


def x_faces(s):
    if "x_face" in s:
        return np.asarray(s["x_face"])
    xf = [-1.0]  # 旧格式：由格中心反推格面
    for xc in np.asarray(s["x"]):
        xf.append(2 * xc - xf[-1])
    return np.array(xf)


def phase(run, s, sp, side):
    """返回 (x 边界, v 中心, v 边界, F[ix, iv], x 标签, v 标签)；对折时 x 为离壁距离 (λ)，v 指向壁面为正。"""
    key = f"f_{sp}"
    if key not in s:
        raise KeyError(f"{rname(run)}: 快照中没有 {key}")
    f = np.asarray(s[key], float)
    v = np.asarray(s[f"v_{sp}"])
    m = mass_of(s, sp)
    if sp == "i":  # 离子用 v_Ti
        v, f = v * np.sqrt(m), f / np.sqrt(m)
    vlab = "v [v_Te]" if sp == "e" else "v [v_Ti]"
    xf = x_faces(s)
    if side == "full":
        X, F, xlab = xf, f, "x [L]"
    else:
        lam = lam_of(run, s)
        n = len(xf) // 2  # 每侧格数
        X = (1.0 - xf[-n - 1:])[::-1] / lam
        right = f[-n:][::-1]
        left = f[:n][:, ::-1]  # 左壁: 指向壁面为 −v，速度格对称，翻转即 v → −v
        F = {"right": right, "left": left, "avg": 0.5 * (left + right)}[side]
        xlab = "distance from wall [λ]"
        vlab = vlab.replace("v [", "v (toward wall +) [")
    dv = v[1] - v[0]
    ve = np.concatenate([v - 0.5 * dv, [v[-1] + 0.5 * dv]])
    return X, v, ve, F, xlab, vlab


def crop(X, F, args):
    """只保留 dmax 以内的格。"""
    if args.side != "full" and args.dmax:
        n = int(np.searchsorted(X, args.dmax)) + 1
        return X[: n + 1], F[:n]
    return X, F


def auto_vlim(v, Fs, floor):
    fmax = max(F.max() for F in Fs)
    hit = np.zeros(len(v), bool)
    for F in Fs:
        hit |= (F > floor * fmax).any(axis=0)
    return 1.15 * np.abs(v[hit]).max() if hit.any() else np.abs(v).max()


def make_norm(Fs, args):
    fmax = max(F.max() for F in Fs)
    if args.lin:
        return Normalize(0, fmax)
    return LogNorm(vmin=args.floor * fmax, vmax=fmax)


def draw(ax, X, ve, F, norm, args, logd_ok=True):
    Z = F.T
    if isinstance(norm, LogNorm):
        Z = np.maximum(Z, norm.vmin)
    if args.logd and args.side != "full" and logd_ok:
        X = X.copy()
        X[0] = 0.5 * X[1]
        ax.set_xscale("log")
    pc = ax.pcolormesh(X, ve, Z, norm=norm, cmap=SEQ, shading="flat", rasterized=True)
    ax.grid(False)
    return pc


def species_in(run, sps):
    """去掉快照中没有的粒子 (旧格式只有电子)。"""
    s = run.snap(0)
    keep = [sp for sp in sps if f"f_{sp}" in s]
    for sp in set(sps) - set(keep):
        print(f"{rname(run)}: 没有 f_{sp}，跳过")
    if not keep:
        sys.exit("没有可画的分布函数")
    return keep


def pick_times(run, tstr, nmax=4):
    ts = run.snap_t
    if tstr:
        out = []
        for t in (float(v) for v in tstr.split(",")):
            k = int(np.argmin(np.abs(ts - t)))
            if abs(ts[k] - t) > 0.01:
                print(f"{rname(run)}: 没有 t={t:g} 的快照，用最近的 t={ts[k]:.4g}")
            out.append(k)
        return out
    return list(np.unique(np.linspace(0, len(ts) - 1, nmax).round().astype(int)))


def save(fig, args, default):
    out = args.out or default
    fig.savefig(out, dpi=args.dpi, bbox_inches="tight")
    print("图已保存:", out)
    if args.show:
        plt.show()


# ======================= 子命令 =======================
def cmd_panel(args):
    run = Run(args.runs[0])
    ks = pick_times(run, args.t)
    sps = species_in(run, args.sp.split(","))
    snaps = [run.snap(k) for k in ks]
    fig, axs = plt.subplots(len(sps), len(ks), figsize=(3.3 * len(ks) + 0.8, 2.8 * len(sps)), squeeze=False,
                            sharey="row", constrained_layout=True)
    for r, sp in enumerate(sps):
        data = []
        for s in snaps:
            X, v, ve, F, xlab, vlab = phase(run, s, sp, args.side)
            X, F = crop(X, F, args)
            data.append((X, F))
        norm = make_norm([F for _, F in data], args)
        vlim = args.vlim or auto_vlim(v, [F for _, F in data], args.floor)
        for c, ((X, F), s) in enumerate(zip(data, snaps)):
            ax = axs[r][c]
            pc = draw(ax, X, ve, F, norm, args)
            ax.set_ylim(-vlim, vlim)
            ax.axhline(0, color=GRAY, lw=0.5)
            ax.set_title(f"f_{sp}, t = {float(s['t']):.3g} τi")
            if r == len(sps) - 1:
                ax.set_xlabel(xlab)
            if c == 0:
                ax.set_ylabel(vlab)
        fig.colorbar(pc, ax=axs[r].tolist(), shrink=0.9).set_label(f"f_{sp}")
    fig.suptitle(f"{rname(run)} ({args.side})", fontsize=10)
    save(fig, args, os.path.join(run.path, f"dist_panel_{'_'.join(sps)}.png"))


def cmd_anim(args):
    run = Run(args.runs[0])
    sps = species_in(run, args.sp.split(","))
    ks = list(range(len(run.snap_files)))
    if args.tmax is not None:
        ks = [k for k in ks if run.snap_t[k] <= args.tmax + 1e-9]
    # 色标与速度范围取全部帧
    frames = {sp: [] for sp in sps}
    for k in ks:
        s = run.snap(k)
        for sp in sps:
            X, v, ve, F, xlab, vlab = phase(run, s, sp, args.side)
            X, F = crop(X, F, args)
            frames[sp].append((X, F, float(s["t"]), v, ve, xlab, vlab))
    fig, axs = plt.subplots(1, len(sps), figsize=(5.0 * len(sps), 3.8), squeeze=False, constrained_layout=True)
    meshes, titles = [], []
    for c, sp in enumerate(sps):
        fr = frames[sp]
        norm = make_norm([F for _, F, *_ in fr], args)
        v, ve, xlab, vlab = fr[0][3], fr[0][4], fr[0][5], fr[0][6]
        vlim = args.vlim or auto_vlim(v, [F for _, F, *_ in fr], args.floor)
        ax = axs[0][c]
        pc = draw(ax, fr[0][0], ve, fr[0][1], norm, args)
        ax.set_ylim(-vlim, vlim)
        ax.set_xlabel(xlab)
        ax.set_ylabel(vlab)
        fig.colorbar(pc, ax=ax).set_label(f"f_{sp}")
        meshes.append((pc, norm))
        titles.append(ax.set_title(""))

    def update(i):
        for c, sp in enumerate(sps):
            X, F, t, *_ = frames[sp][i]
            pc, norm = meshes[c]
            Z = F.T if not isinstance(norm, LogNorm) else np.maximum(F.T, norm.vmin)
            pc.set_array(Z.ravel())
            titles[c].set_text(f"{rname(run)}: f_{sp}, t = {t:.3g} τi")
        return [m for m, _ in meshes] + titles

    update(0)
    ani = FuncAnimation(fig, update, frames=len(ks), interval=1000 / args.fps)
    out = args.out or os.path.join(run.path, f"dist_anim_{'_'.join(sps)}.gif")
    ani.save(out, writer=PillowWriter(fps=args.fps), dpi=args.dpi // 2 + 20)
    print("动画已保存:", out, f"({len(ks)} 帧)")


def cmd_slice(args):
    runs = [Run(r) for r in args.runs]
    ds = [float(v) for v in (args.d or "0.25,10,100").split(",")]
    sp = args.sp.split(",")[0]
    multi = len(runs) > 1
    fig, axs = plt.subplots(1, len(ds), figsize=(3.8 * len(ds), 3.3), squeeze=False, sharey=True, constrained_layout=True)
    vmax_all = 0.0
    for ri, run in enumerate(runs):
        ks = pick_times(run, args.t, nmax=3)
        for ti, k in enumerate(ks):
            s = run.snap(k)
            if f"f_{sp}" not in s:
                print(f"{rname(run)}: 没有 f_{sp} (旧格式只有电子)，跳过")
                break
            X, v, ve, F, xlab, vlab = phase(run, s, sp, args.side if args.side != "full" else "right")
            dc = 0.5 * (X[1:] + X[:-1])
            if multi:
                color, ls = CAT[ri % len(CAT)], DASH[ti % len(DASH)]
            else:
                color, ls = BLUE[int(round(ti * (len(BLUE) - 1) / max(len(ks) - 1, 1)))], "-"
            for di, d in enumerate(ds):
                i = int(np.argmin(np.abs(dc - d)))
                ax = axs[0][di]
                fv = F[i]
                lab = (f"{rname(run)} " if multi else "") + f"t={float(s['t']):.3g}"
                if multi:  # 不同算例网格可能不同：标出实际格中心距离
                    lab += f", d={dc[i]:.3g}λ"
                ax.plot(v, fv, color=color, ls=ls, label=lab)
                if args.maxw:
                    n, u, T, _ = moments(fv, v, 1.0)  # v 已归一到本粒子热速度, 质量取 1
                    if n > 0 and T > 0:
                        ax.plot(v, n / np.sqrt(2 * np.pi * T) * np.exp(-(v - u) ** 2 / (2 * T)),
                                color=GRAY, lw=0.9, ls=(0, (3, 2)))
                ax.set_title(f"d ≈ {d:g} λ" if multi else f"d = {dc[i]:.3g} λ")
                vmax_all = max(vmax_all, fv.max())
                hit = np.abs(v[fv > args.floor * fv.max()]) if fv.max() > 0 else np.array([0.0])
                ax._vl = max(getattr(ax, "_vl", 0.0), 1.15 * hit.max())
    for ax in axs[0]:
        ax.set_xlabel(vlab)
        lim = args.vlim or ax._vl
        ax.set_xlim(-lim, lim)
        if not args.lin:
            ax.set_yscale("log")
            ax.set_ylim(args.floor * vmax_all, 2 * vmax_all)
    axs[0][0].set_ylabel(f"f_{sp}")
    h, lab = axs[0][0].get_legend_handles_labels()
    if args.maxw:
        h.append(plt.Line2D([], [], color=GRAY, lw=0.9, ls=(0, (3, 2))))
        lab.append("local Maxwellian")
    axs[0][0].legend(h, lab, fontsize=7)
    save(fig, args, os.path.join(runs[0].path, f"dist_slice_{sp}.png"))


def cmd_vt(args):
    run = Run(args.runs[0])
    sp = args.sp.split(",")[0]
    d = float((args.d or "0.25").split(",")[0])
    rows, ts = [], []
    for k in range(len(run.snap_files)):
        if args.tmax is not None and run.snap_t[k] > args.tmax + 1e-9:
            continue
        s = run.snap(k)
        X, v, ve, F, xlab, vlab = phase(run, s, sp, args.side if args.side != "full" else "right")
        dc = 0.5 * (X[1:] + X[:-1])
        i = int(np.argmin(np.abs(dc - d)))
        rows.append(F[i])
        ts.append(float(s["t"]))
    Z, t = np.array(rows), np.array(ts)
    if len(t) < 2:
        sys.exit("快照少于 2 个")
    te = np.concatenate([[t[0] - 0.5 * (t[1] - t[0])], 0.5 * (t[1:] + t[:-1]), [t[-1] + 0.5 * (t[-1] - t[-2])]])
    norm = make_norm([Z], args)
    vlim = args.vlim or auto_vlim(v, [Z], args.floor)
    fig, ax = plt.subplots(figsize=(6.4, 4.0), constrained_layout=True)
    pc = ax.pcolormesh(ve, te, np.maximum(Z, norm.vmin) if isinstance(norm, LogNorm) else Z,
                       norm=norm, cmap=SEQ, shading="flat", rasterized=True)
    fig.colorbar(pc, ax=ax).set_label(f"f_{sp}")
    ax.set_xlim(-vlim, vlim)
    ax.set_xlabel(vlab)
    ax.set_ylabel("t [τi]")
    ax.grid(False)
    ax.set_title(f"{rname(run)}: f_{sp} at d = {dc[i]:.3g} λ ({args.side})")
    save(fig, args, os.path.join(run.path, f"dist_vt_{sp}_d{dc[i]:.3g}.png"))


def cmd_diff(args):
    if len(args.runs) != 2:
        sys.exit("diff 需要两个算例: RUN RUN0 (画 RUN − RUN0)")
    run, ref = Run(args.runs[0]), Run(args.runs[1])
    sps = args.sp.split(",")
    k = pick_times(run, args.t or "0.5", 1)[0]
    s = run.snap(k)
    s0 = ref.snap(t=float(s["t"]))
    if abs(float(s0["t"]) - float(s["t"])) > 0.01:
        print(f"{rname(ref)}: 最近的快照 t={float(s0['t']):.4g}，与 {float(s['t']):.4g} 不同")
    fig, axs = plt.subplots(1, len(sps), figsize=(5.0 * len(sps), 3.8), squeeze=False, constrained_layout=True)
    for c, sp in enumerate(sps):
        X, v, ve, F, xlab, vlab = phase(run, s, sp, args.side)
        X0, _, _, F0, _, _ = phase(ref, s0, sp, args.side)
        if F.shape != F0.shape:
            sys.exit(f"网格不同: {F.shape} vs {F0.shape}")
        X, D = crop(X, F - F0, args)
        _, Fc = crop(X0, F0, args)
        lim = np.abs(D).max() or 1.0
        ax = axs[0][c]
        pc = draw(ax, X, ve, D, TwoSlopeNorm(0.0, -lim, lim), args)
        pc.set_cmap(DIV)
        vlim = args.vlim or auto_vlim(v, [Fc], args.floor)
        ax.set_ylim(-vlim, vlim)
        ax.set_xlabel(xlab)
        ax.set_ylabel(vlab)
        ax.set_title(f"f_{sp}: {rname(run)} − {rname(ref)}, t = {float(s['t']):.3g} τi")
        fig.colorbar(pc, ax=ax).set_label(f"Δf_{sp} (max |f0| = {np.abs(Fc).max():.3g})")
    save(fig, args, os.path.join(run.path, f"dist_diff_{rname(ref)}.png"))


def cell_avg(X, F, a, b):
    """折叠坐标下 [a, b] 区间内各格 f 的格宽加权平均 (把 REF 限制到 TEST 的格上)。"""
    lo, hi = np.maximum(X[:-1], a), np.minimum(X[1:], b)
    w = np.clip(hi - lo, 0.0, None)
    return (w[:, None] * F).sum(axis=0) / w.sum()


def cmd_wall(args):
    """TEST 与 REF 在给定离壁距离处 f(v) 的逐时刻对比；REF 按 TEST 的格做格宽加权平均。"""
    if len(args.runs) != 2:
        sys.exit("wall 需要两个算例: TEST REF")
    run, ref = Run(args.runs[0]), Run(args.runs[1])
    sp = args.sp.split(",")[0]
    d = float((args.d or "0.25").split(",")[0])
    side = args.side if args.side != "full" else "avg"
    ks = pick_times(run, args.t or "0.2,0.4,1,2,3,4")
    ncol = min(3, len(ks))
    nrow = (len(ks) + ncol - 1) // ncol
    fig = plt.figure(figsize=(3.9 * ncol, 3.9 * nrow), constrained_layout=True)
    gs = fig.add_gridspec(2 * nrow, ncol, height_ratios=[3, 1.3] * nrow)
    cT, cR = CAT[0], CAT[1]
    for n, k in enumerate(ks):
        r, c = divmod(n, ncol)
        s = run.snap(k)
        t = float(s["t"])
        k0 = int(np.argmin(np.abs(ref.snap_t - t)))
        s0 = ref.snap(k0)
        X, v, _, F, xlab, vlab = phase(run, s, sp, side)
        X0, v0, _, F0, _, _ = phase(ref, s0, sp, side)
        if len(v) != len(v0):
            sys.exit(f"速度网格不同: {len(v)} vs {len(v0)}")
        i = int(np.argmin(np.abs(0.5 * (X[1:] + X[:-1]) - d)))
        a, b = X[i], X[i + 1]
        fT, fR = F[i], cell_avg(X0, F0, a, b)
        mT, mR = moments(fT, v, 1.0), moments(fR, v, 1.0)
        ax = fig.add_subplot(gs[2 * r, c])
        ax.plot(v, fR, color=cR, lw=1.8, label=f"{rname(ref)} (avg over cell)")
        ax.plot(v, fT, color=cT, lw=1.1, ls=(0, (4, 2)), label=rname(run))
        fmax = max(fT.max(), fR.max())
        hit = np.abs(v[(fR > args.floor * fmax) | (fT > args.floor * fmax)])
        vl = args.vlim or 1.15 * hit.max()
        ax.set_xlim(-vl, vl)
        if not args.lin:
            ax.set_yscale("log")
            ax.set_ylim(args.floor * fmax, 2 * fmax)
        dn = (mT[0] / mR[0] - 1) * 100
        dT = (mT[2] / mR[2] - 1) * 100
        ax.set_title(f"t = {t:.3g} τi   Δn {dn:+.2f}%  ΔT {dT:+.2f}%", fontsize=9)
        if c == 0:
            ax.set_ylabel(f"f_{sp}")
        if n == 0:
            ax.legend(fontsize=7, loc="lower center")
        ax.tick_params(labelbottom=False)
        bx = fig.add_subplot(gs[2 * r + 1, c], sharex=ax)
        bx.axhline(0, color=GRAY, lw=0.6)
        bx.plot(v, (fT - fR) / fR.max() * 100, color=cT, lw=1.0)
        bx.set_xlabel(vlab)
        if c == 0:
            bx.set_ylabel("Δf / max f_REF [%]")
    fig.suptitle(f"f_{sp} at {d:g} λ from wall (cell {a:.3g}–{b:.3g} λ, side={side}): "
                 f"{rname(run)} vs {rname(ref)}", fontsize=10)
    save(fig, args, os.path.join(run.path, f"dist_wall_{sp}_{rname(ref)}.png"))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["panel", "anim", "slice", "vt", "diff", "wall"])
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--t", help="时刻 (τi)，逗号分隔")
    ap.add_argument("--sp", default="e,i", help="粒子: e、i 或 e,i")
    ap.add_argument("--d", help="离壁距离 (λ)，逗号分隔")
    ap.add_argument("--side", default="avg", choices=["avg", "left", "right", "full"])
    ap.add_argument("--dmax", type=float)
    ap.add_argument("--logd", action="store_true")
    ap.add_argument("--lin", action="store_true", help="线性色标 / 线性纵轴")
    ap.add_argument("--floor", type=float, default=1e-6, help="对数色标下限 (相对最大值)")
    ap.add_argument("--vlim", type=float)
    ap.add_argument("--tmax", type=float)
    ap.add_argument("--maxw", action="store_true", help="slice: 叠加局部麦克斯韦分布")
    ap.add_argument("--fps", type=float, default=4)
    ap.add_argument("--out")
    ap.add_argument("--dpi", type=int, default=150)
    ap.add_argument("--show", action="store_true")
    args = ap.parse_args()
    {"panel": cmd_panel, "anim": cmd_anim, "slice": cmd_slice, "vt": cmd_vt, "diff": cmd_diff, "wall": cmd_wall}[args.cmd](args)


if __name__ == "__main__":
    main()
