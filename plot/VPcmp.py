"""两个算例逐点对比：把参考算例 (细网格) 守恒地限制到被检算例 (粗 / 非均匀网格) 的格上，再逐格比较。

用法:
  python plot/VPcmp.py TEST REF [--t 0.5,1,2,4] [--d 0.25,4,20,100] [--side avg] [--out 前缀]
例:
  python plot/VPcmp.py runs/Case84 runs/Case83          # refine 280 格 对 均匀 4000 格

要求 TEST 的每个格面都是 REF 的格面 (refine / uniform 嵌套网格满足)，两者速度网格相同、Δt 相同。
限制方式 (格平均，守恒):
  n、ρ、φ、E、f 直接按格宽加权平均；u、T、q 由限制后的 ∫f v^k dv (k = 0…3) 重组，
  即与"对 REF 的 f 先做格平均再求矩"完全相同。E、φ 在 TEST 中是格中心值，与格平均差 O(Δx² E'')。
输出 (存在 TEST 目录下，前缀 cmp_<REF>):
  _summary.txt   各时刻、各区域 (贴壁细格块 / 其后 4 个粗格 / 至 200λ / 体区) 的最大偏差
  _prof.png      各时刻剖面叠画 (REF 原始、REF 限制、TEST) 与逐格相对偏差
  _err_t.png     各量最大偏差随时间变化 (分区域)
  _f.png         f 的逐格偏差 (相空间) 与 f 的 L1 相对误差随时间变化
  _diag.png      壁面通量、能流、壁面 n/T 的相对偏差随时间变化
偏差定义: n、T 为逐点相对偏差 |a−b|/b；有符号量 (u、q、E、φ、ρ) 为 |a−b| / 该区域内 max|b|。
旧格式算例 (仅电子 f、每 OutF 一个快照) 也可对比，只比较电子。
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vpio import Run  # noqa: E402

import matplotlib  # noqa: E402

if "--show" not in sys.argv:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import TwoSlopeNorm  # noqa: E402

from VPplot import BLUE, CAT, DIV, LABEL, SIGNED, VECTOR, Source  # noqa: E402

GRAY = "#8a8984"
POS = ["n_e", "n_i", "T_e", "T_i"]
SGN = ["u_e", "u_i", "q_e", "q_i", "phi", "E", "rho"]


def faces(src):
    return np.concatenate([[src.x[0] - 0.5 * src.dx[0]], src.x + 0.5 * src.dx])


class Restrict:
    """REF 网格 → TEST 网格的守恒限制。"""

    def __init__(self, fine, coarse):
        xf, xc = faces(fine), faces(coarse)
        idx = np.searchsorted(xf, xc - 1e-9 * fine.dx.min())
        idx = np.clip(idx, 0, len(xf) - 1)
        bad = np.abs(xf[idx] - xc) > 1e-6 * fine.dx.min()
        if bad.any():
            sys.exit(f"{coarse.name} 的格面不全是 {fine.name} 的格面 (第一个不符: x = {xc[bad][0]:.9g})，无法逐格限制")
        self.start = idx[:-1]
        self.wf = fine.dx
        self.wc = np.add.reduceat(fine.dx, self.start)
        self.ratio = np.diff(idx)

    def __call__(self, q):
        q = np.asarray(q, float)
        w = self.wf if q.ndim == 1 else self.wf[:, None]
        wc = self.wc if q.ndim == 1 else self.wc[:, None]
        return np.add.reduceat(q * w, self.start, axis=0) / wc


def restrict_profiles(R, d, mass):
    """限制一组剖面；u、T、q 由矩重组。"""
    out = {}
    for k in ("phi", "E", "rho"):
        if k in d:
            out[k] = R(d[k])
    for sp, m in mass.items():
        if f"n_{sp}" not in d:
            continue
        n, u, T, q = (np.asarray(d[f"{a}_{sp}"], float) for a in "nuTq")
        M1 = n * u
        M2 = n * T / m + n * u**2
        M3 = 2 * q / m + 3 * u * n * T / m + n * u**3
        nc, M1c, M2c, M3c = R(n), R(M1), R(M2), R(M3)
        uc = np.where(nc > 0, M1c / np.where(nc > 0, nc, 1), 0)
        Tm = np.where(nc > 0, M2c / np.where(nc > 0, nc, 1) - uc**2, 0)
        out[f"n_{sp}"], out[f"u_{sp}"], out[f"T_{sp}"] = nc, uc, m * Tm
        out[f"q_{sp}"] = 0.5 * m * (M3c - 3 * uc * nc * Tm - nc * uc**3)
    return out


def regions(src):
    """离壁距离 (λ) 的分区：贴壁细格块、其后 4 个粗格、至 200λ、体区。"""
    itf = src.interfaces()
    w = itf[0] if len(itf) else 8.0
    n = len(src.x) // 2
    h = src.dx[-n:][::-1] / src.lam
    coarse = h[np.searchsorted(np.cumsum(h), w + 1e-9)] if w < np.sum(h) else h[-1]
    b2 = w + 4 * coarse
    return [(f"0–{w:g}λ", 0, w), (f"{w:g}–{b2:g}λ", w, b2), (f"{b2:g}–200λ", b2, 200), ("200λ–center", 200, 1e9)]


def errors(a, b, name, dist, regs):
    """各区域的最大偏差。"""
    out = []
    for _, lo, hi in regs:
        m = (dist >= lo) & (dist < hi)
        if not m.any():
            out.append(np.nan)
            continue
        if name in POS:
            e = np.abs(a[m] - b[m]) / np.maximum(np.abs(b[m]), 1e-300)
        else:
            e = np.abs(a[m] - b[m]) / max(np.abs(b[m]).max(), 1e-300)
        out.append(e.max())
    return out


def match_times(test, ref):
    """按 ts 配对 (同一 Δt)；否则按最近时刻。"""
    pairs = []
    for k, t in enumerate(test.t):
        kr = int(np.argmin(np.abs(ref.t - t)))
        if abs(ref.t[kr] - t) < 1e-6 + 1e-3 * max(t, 1e-3):
            pairs.append((k, kr))
    return pairs


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("test")
    ap.add_argument("ref")
    ap.add_argument("--t", default="0.5,1,2,4", help="剖面图与 f 图的时刻 (τi)")
    ap.add_argument("--vars", default="E,phi,n_e,T_e,u_i,T_i", help="剖面图的变量")
    ap.add_argument("--side", default="avg", choices=["avg", "left", "right"])
    ap.add_argument("--dmax", type=float, default=1000.0)
    ap.add_argument("--out", help="输出文件前缀 (缺省 TEST/cmp_<REF>)")
    ap.add_argument("--dpi", type=int, default=150)
    ap.add_argument("--show", action="store_true")
    args = ap.parse_args()

    test, ref = Source(args.test), Source(args.ref)
    R = Restrict(ref, test)
    regs = regions(test)
    dist = test.fold(test.x, "x", args.side)[0]
    pre = args.out or os.path.join(test.path, f"cmp_{ref.name}")
    mass = {"e": 1.0, "i": 1836.0}
    pairs = match_times(test, ref)
    if not pairs:
        sys.exit("两个算例没有相同的输出时刻")
    print(f"{test.name} ({len(test.x)} 格) 对 {ref.name} ({len(ref.x)} 格，限制比 {R.ratio.min()}–{R.ratio.max()})，"
          f"共同时刻 {len(pairs)} 个，数据来源: {'剖面 prof.vpd' if test.has_prof and ref.has_prof else '快照的矩'}")

    # ---------- 剖面: 各量、各区域的最大偏差随时间 ----------
    names = None
    hist = {}
    times = []
    lines = []
    for k, kr in pairs:
        dt_, dr = test.raw(k), ref.raw(kr)
        rr = restrict_profiles(R, dr, mass)
        if names is None:
            names = [n for n in POS + SGN if n in dt_ and n in rr]
        times.append(test.t[k])
        for n in names:
            a = test.fold(np.asarray(dt_[n], float), n, args.side)[1]
            b = test.fold(rr[n], n, args.side)[1]
            hist.setdefault(n, []).append(errors(a, b, n, dist, regs))
    times = np.array(times)
    hist = {n: np.array(v) for n, v in hist.items()}

    head = f"{'t':>6} {'量':>5} " + " ".join(f"{r[0]:>14}" for r in regs)
    tshow = [float(v) for v in args.t.split(",")]
    sel = sorted({int(np.argmin(np.abs(times - t))) for t in tshow})
    lines.append(f"{test.name} 对 {ref.name} (REF 限制到 TEST 网格)，各区域最大偏差 [n、T: 逐点相对；其余: 相对区域内 max|REF|]")
    lines.append(head)
    for j in sel:
        for n in names:
            lines.append(f"{times[j]:6.3g} {n:>5} " + " ".join(f"{100 * e:13.3g}%" for e in hist[n][j]))
        lines.append("")
    lines.append("全时段最大 (t ≥ 0.1 τi):")
    late = times >= 0.1
    for n in names:
        mx = np.nanmax(hist[n][late], axis=0) if late.any() else hist[n].max(axis=0)
        lines.append(f"{'':6} {n:>5} " + " ".join(f"{100 * e:13.3g}%" for e in mx))

    # ---------- 剖面图 ----------
    pv = [v for v in args.vars.split(",") if v in names]
    if pv:
        fig, axs = plt.subplots(2, len(pv), figsize=(3.6 * len(pv), 5.2), squeeze=False, sharex=True,
                                gridspec_kw={"height_ratios": [2, 1]})
        for ti, t in enumerate(tshow):
            j = int(np.argmin(np.abs(times - t)))
            k, kr = pairs[j]
            dt_, dr = test.raw(k), ref.raw(kr)
            rr = restrict_profiles(R, dr, mass)
            col = BLUE[int(round(ti * (len(BLUE) - 1) / max(len(tshow) - 1, 1)))]
            for c, n in enumerate(pv):
                xr, br = ref.fold(np.asarray(dr[n], float), n, args.side)
                xt, a = test.fold(np.asarray(dt_[n], float), n, args.side)
                _, b = test.fold(rr[n], n, args.side)
                ax = axs[0][c]
                ax.plot(xr, br, color=col, lw=1.0, alpha=0.6)
                ax.plot(xt, a, color=col, lw=0, marker="o", ms=3.5, label=f"t={times[j]:.3g}")
                if n in POS:
                    rel = (a - b) / np.maximum(np.abs(b), 1e-300)
                else:
                    rel = (a - b) / max(np.abs(b).max(), 1e-300)
                axs[1][c].plot(xt, 100 * rel, color=col, marker="o", ms=2.5, lw=1.0)
        for c, n in enumerate(pv):
            axs[0][c].set_ylabel(LABEL.get(n, n))
            axs[1][c].set_ylabel("Δ [%]" + ("" if n in POS else " of max|REF|"))
            axs[1][c].axhline(0, color=GRAY, lw=0.6)
            axs[1][c].set_xlabel("distance from wall [λ]")
            for ax in (axs[0][c], axs[1][c]):
                ax.set_xscale("log")
                ax.set_xlim(0.15, args.dmax)
                for d in test.interfaces():
                    ax.axvline(d, color=GRAY, lw=0.8, ls=(0, (2, 2)), zorder=0)
            if n in SIGNED:
                axs[0][c].axhline(0, color=GRAY, lw=0.6)
        h = [plt.Line2D([], [], color=GRAY, lw=1.0), plt.Line2D([], [], color=GRAY, lw=0, marker="o", ms=3.5)]
        h += [plt.Line2D([], [], color=BLUE[int(round(i * (len(BLUE) - 1) / max(len(tshow) - 1, 1)))], lw=2)
              for i in range(len(tshow))]
        axs[0][0].legend(h, [f"{ref.name} (raw)", test.name] + [f"t={t:g}" for t in tshow], fontsize=7)
        fig.suptitle(f"{test.name} vs {ref.name}: Δ = TEST − REF restricted to TEST cells ({args.side})", fontsize=10)
        fig.tight_layout()
        fig.savefig(pre + "_prof.png", dpi=args.dpi, bbox_inches="tight")
        print("图已保存:", pre + "_prof.png")
        plt.close(fig)

    # ---------- 偏差随时间 ----------
    fig, axs = plt.subplots(1, len(regs), figsize=(3.6 * len(regs), 3.2), squeeze=False, sharey=True)
    for r, (rn, _, _) in enumerate(regs):
        ax = axs[0][r]
        for c, n in enumerate(names):
            ax.plot(times, 100 * hist[n][:, r], color=CAT[c % len(CAT)], ls="-" if c < len(CAT) else "--", label=n)
        ax.set_yscale("log")
        ax.set_title(rn)
        ax.set_xlabel("t [τi]")
    axs[0][0].set_ylabel("max deviation [%]")
    axs[0][-1].legend(fontsize=7, ncol=2)
    fig.suptitle(f"{test.name} vs {ref.name}", fontsize=10)
    fig.tight_layout()
    fig.savefig(pre + "_err_t.png", dpi=args.dpi, bbox_inches="tight")
    print("图已保存:", pre + "_err_t.png")
    plt.close(fig)

    # ---------- 分布函数 ----------
    rt, rr_ = Run(args.test), Run(args.ref)
    if rt.snap_files and rr_.snap_files:
        st, sr = rt.snap_t, rr_.snap_t
        fpairs = [(k, int(np.argmin(np.abs(sr - t)))) for k, t in enumerate(st)]
        fpairs = [(k, kr) for k, kr in fpairs if abs(sr[kr] - st[k]) < 1e-6 + 1e-3 * max(st[k], 1e-3)]
        sps = [sp for sp in ("e", "i") if f"f_{sp}" in rt.snap(0) and f"f_{sp}" in rr_.snap(0)]
        fl1 = {sp: [] for sp in sps}
        ft = []
        figsel = sorted({int(np.argmin(np.abs(np.array([st[k] for k, _ in fpairs]) - t))) for t in tshow})
        fig, axs = plt.subplots(len(sps), len(figsel), figsize=(3.4 * len(figsel) + 0.8, 2.9 * len(sps)),
                                squeeze=False, constrained_layout=True)
        dvc = {}
        for jj, (k, kr) in enumerate(fpairs):
            a_s, b_s = rt.snap(k), rr_.snap(kr)
            ft.append(st[k])
            for r_, sp in enumerate(sps):
                va, vb = np.asarray(a_s[f"v_{sp}"]), np.asarray(b_s[f"v_{sp}"])
                if va.shape != vb.shape or np.abs(va - vb).max() > 1e-9:
                    sys.exit("速度网格不同，无法逐点比较 f")
                fa = np.asarray(a_s[f"f_{sp}"], float)
                fb = R(np.asarray(b_s[f"f_{sp}"], float))
                w = test.dx[:, None]
                dvc[sp] = va[1] - va[0]
                e = []
                for _, lo, hi in regs:
                    dd = (1 - np.abs(test.x)) / test.lam
                    m = (dd >= lo) & (dd < hi)
                    e.append(np.sum(np.abs(fa[m] - fb[m]) * w[m]) / max(np.sum(np.abs(fb[m]) * w[m]), 1e-300))
                fl1[sp].append(e)
                if jj in figsel:
                    c = figsel.index(jj)
                    n = len(test.x) // 2
                    D = 0.5 * ((fa - fb)[-n:][::-1] + (fa - fb)[:n][:, ::-1])  # 两壁对折平均，v 指向壁面为正
                    F0 = 0.5 * (fb[-n:][::-1] + fb[:n][:, ::-1])
                    X = (1.0 - faces(test)[-n - 1:])[::-1] / test.lam
                    X[0] = 0.5 * X[1]
                    vv = va * (np.sqrt(mass[sp]) if sp == "i" else 1.0)
                    dv = vv[1] - vv[0]
                    ve = np.concatenate([vv - 0.5 * dv, [vv[-1] + 0.5 * dv]])
                    lim = np.abs(D).max() or 1.0
                    ax = axs[r_][c]
                    pc = ax.pcolormesh(X, ve, D.T, cmap=DIV, norm=TwoSlopeNorm(0, -lim, lim), shading="flat", rasterized=True)
                    hit = np.abs(vv[(F0 > 1e-6 * F0.max()).any(axis=0)])
                    vl = 1.15 * hit.max() if hit.size else np.abs(vv).max()
                    ax.set_ylim(-vl, vl)
                    ax.set_xscale("log")
                    ax.set_xlim(right=args.dmax)
                    ax.grid(False)
                    ax.set_title(f"Δf_{sp}, t={st[k]:.3g} (max f={F0.max():.3g})", fontsize=8)
                    if c == 0:
                        ax.set_ylabel("v (toward wall +) " + ("[v_Te]" if sp == "e" else "[v_Ti]"))
                    if r_ == len(sps) - 1:
                        ax.set_xlabel("distance from wall [λ]")
                    fig.colorbar(pc, ax=ax)
        fig.suptitle(f"f: {rt.path.rstrip('/').split('/')[-1]} − {rr_.path.rstrip('/').split('/')[-1]} (restricted)", fontsize=10)
        fig.savefig(pre + "_f.png", dpi=args.dpi, bbox_inches="tight")
        print("图已保存:", pre + "_f.png")
        plt.close(fig)
        ft = np.array(ft)
        lines.append("")
        lines.append("分布函数 L1 相对误差 Σ|Δf|ΔxΔv / Σ|f_REF|ΔxΔv (两壁合计):")
        lines.append(f"{'t':>6} {'粒子':>4} " + " ".join(f"{r[0]:>14}" for r in regs))
        for sp in sps:
            arr = np.array(fl1[sp])
            for j, t in enumerate(ft):
                lines.append(f"{t:6.3g} {sp:>4} " + " ".join(f"{100 * e:13.3g}%" for e in arr[j]))
        fig, axs = plt.subplots(1, len(sps), figsize=(4.2 * len(sps), 3.2), squeeze=False, sharey=True)
        for r_, sp in enumerate(sps):
            arr = np.array(fl1[sp])
            for c, (rn, _, _) in enumerate(regs):
                axs[0][r_].plot(ft, 100 * arr[:, c], color=CAT[c], marker="o", ms=3, label=rn)
            axs[0][r_].set_yscale("log")
            axs[0][r_].set_title(f"f_{sp}: L1 relative error")
            axs[0][r_].set_xlabel("t [τi]")
        axs[0][0].set_ylabel("[%]")
        axs[0][0].legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(pre + "_fL1.png", dpi=args.dpi, bbox_inches="tight")
        print("图已保存:", pre + "_fL1.png")
        plt.close(fig)

    # ---------- 壁面诊断 ----------
    try:
        da, db = rt.diag, rr_.diag
        cols = [c for c in ("flux.H+", "flux.e-", "Enflux.H+", "Enflux.e-", "nw.H+", "Tw.H+", "nw.e-", "Tw.e-") if c in da and c in db]
        fig, axs = plt.subplots(2, 1, figsize=(6.4, 5.0), sharex=True)
        lines.append("")
        lines.append("壁面量相对偏差 (TEST / REF − 1):")
        lines.append(f"{'t':>6} " + " ".join(f"{c:>10}" for c in cols))
        tq = [t for t in (0.25, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4) if t <= min(da["t"][-1], db["t"][-1]) + 1e-6]
        rel = {}
        for c in cols:
            bi = np.interp(da["t"], db["t"], db[c])
            rel[c] = np.where(np.abs(bi) > 0, da[c] / np.where(np.abs(bi) > 0, bi, 1) - 1, np.nan)
        for t in tq:
            i = int(np.argmin(np.abs(da["t"] - t)))
            lines.append(f"{t:6.3g} " + " ".join(f"{100 * rel[c][i]:9.3g}%" for c in cols))
        for c_, c in enumerate(cols):
            ax = axs[0] if c.startswith(("flux", "Enflux")) else axs[1]
            m = da["t"] > 0.05
            ax.plot(da["t"][m], 100 * rel[c][m], color=CAT[c_ % len(CAT)], label=c, lw=1.2)
        for ax in axs:
            ax.axhline(0, color=GRAY, lw=0.6)
            ax.set_ylabel("TEST/REF − 1 [%]")
            ax.legend(fontsize=7, ncol=2)
        axs[1].set_xlabel("t [τi]")
        fig.suptitle(f"wall quantities: {test.name} vs {ref.name}", fontsize=10)
        fig.tight_layout()
        fig.savefig(pre + "_diag.png", dpi=args.dpi, bbox_inches="tight")
        print("图已保存:", pre + "_diag.png")
        plt.close(fig)
    except FileNotFoundError:
        pass

    txt = "\n".join(lines)
    open(pre + "_summary.txt", "w").write(txt + "\n")
    print(txt)
    print("汇总已保存:", pre + "_summary.txt")


if __name__ == "__main__":
    main()
