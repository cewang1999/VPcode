"""VPcode 画图 (读 out/*.vpd，旧 txt 算例自动兼容；数据读取见 vpio.py)。

子命令 (RUN 为算例目录，可给多个比较):
  prof RUN... --t 0.3,0.5 --vars E,phi,n_e,T_e   剖面: 离壁距离 (λ) 上的 E、φ、矩等，多时刻 / 多算例叠画
  xt   RUN    --var E [--ref RUN0]                时空图: 离壁距离 × 时间 (可减去对照算例)
  f    RUN    --t 0.5 --sp e [--log]              相空间 f(x, v) (壁面附近)
  fv   RUN... --t 0.5 --sp e --d 0.25,10,100      给定离壁距离处的 f(v)
  diag RUN... --cols ElectricEnergy,Tw.e- [--rel] runtime_diags.csv 时间序列 (--rel 为相对第一个算例)
通用选项: --side avg|left|right|full (缺省 avg: 两壁平均，矢量取指向壁面为正)，
          --dmax 离壁距离上限 (λ)，--logd 距离用对数轴，--out 图文件名，--show 弹窗显示
剖面变量: E phi rho dn(=n_i-n_e) n_e n_i u_e u_i T_e T_i q_e q_i lamD(=局部 λ_D/λ) dxlam(=Δx/λ_D,loc) EE
例: python plot/VPplot.py prof runs/Case78 runs/Case80 --t 0.46,0.5 --vars E,dn,T_e,dxlam --dmax 100 --logd
    python plot/VPplot.py xt runs/Case78 --var E --dmax 100 --logd --tmax 0.7
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vpio import Run, dv_of, moments, v_edges  # noqa: E402

import matplotlib  # noqa: E402

if "--show" not in sys.argv:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, LogNorm, TwoSlopeNorm  # noqa: E402

# ---- 配色 (分类色按固定顺序；有序的多时刻用单色蓝阶；有符号量用蓝-灰-红发散色) ----
CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
BLUE = ["#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]
SEQ = LinearSegmentedColormap.from_list("vp_seq", ["#f4f8fd", "#cde2fb", "#86b6ef", "#3987e5", "#256abf", "#184f95", "#0d366b"])
DIV = LinearSegmentedColormap.from_list("vp_div", ["#0d366b", "#2a78d6", "#9ec5f4", "#f0efec", "#f2a7a6", "#e34948", "#8f2322"])
DASH = ["-", "--", ":", "-."]
VECTOR = {"E", "u_e", "u_i", "q_e", "q_i"}  # 两壁对折时取指向壁面为正的分量
SIGNED = VECTOR | {"rho", "dn"}
LABEL = {
    "E": "E (toward wall +) [T/(eL)]", "phi": "φ [T/e]", "rho": "ρ [e n0]", "dn": "n_i − n_e [n0]",
    "n_e": "n_e [n0]", "n_i": "n_i [n0]", "u_e": "u_e (toward wall +) [v_Te]", "u_i": "u_i (toward wall +) [v_Te]",
    "T_e": "T_e [T0]", "T_i": "T_i [T0]", "q_e": "q_e [n0 T0 v_Te]", "q_i": "q_i [n0 T0 v_Te]",
    "lamD": "local λ_D / λ", "dxlam": "Δx / λ_D,loc",
}

plt.rcParams.update({
    "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb", "axes.edgecolor": "#8a8984",
    "axes.labelcolor": "#0b0b0b", "xtick.color": "#52514e", "ytick.color": "#52514e",
    "axes.grid": True, "grid.color": "#e4e3df", "grid.linewidth": 0.6, "lines.linewidth": 1.6,
    "legend.frameon": False, "font.size": 9, "axes.titlesize": 9, "axes.spines.top": False, "axes.spines.right": False,
})


# ======================= 数据 =======================
class Source:
    """一个算例的剖面数据：优先 out/prof.vpd，其次 VPD 快照 (由 f 求矩)，再次旧 txt (仅电子)。"""

    def __init__(self, path):
        self.path = path.rstrip("/")
        self.name = os.path.basename(self.path)
        self.run = Run(path)
        self.has_prof = os.path.exists(os.path.join(self.run.out, "prof.vpd"))
        if self.has_prof:
            p = self.run.prof
            self.x, self.dx, self.t = np.asarray(p["x"]), np.asarray(p["dx"]), np.asarray(p["t"])
            self.lam = p.meta.get("lambda", 1e-3)
        else:
            if not self.run.snap_files:
                raise FileNotFoundError(f"{path}: 没有 out/prof.vpd，也没有分布函数快照")
            s = self.run.snap(0)
            self.x = np.asarray(s["x"])
            if "dx" in s:
                self.dx = np.asarray(s["dx"])
            else:  # 旧格式只有格中心：由左壁 x = −1 逐格反推格面
                xf = [-1.0]
                for xc in self.x:
                    xf.append(2 * xc - xf[-1])
                self.dx = np.diff(xf)
            self.t = self.run.snap_t
            self.lam = self.run.case.get("lambda", 1e-3)

    def nearest(self, t):
        k = int(np.argmin(np.abs(self.t - t)))
        if abs(self.t[k] - t) > 0.01:
            print(f"{self.name}: 没有 t={t:g} 的数据，用最近的 t={self.t[k]:.4g}")
        return k

    def raw(self, k):
        """第 k 个时刻的全部剖面 (dict)。"""
        if self.has_prof:
            p = self.run.prof
            d = {n: np.asarray(p[n][k]) for n in p.fields if n not in ("x", "x_face", "dx", "v_e", "v_i")}
        else:
            s = self.run.snap(k)
            d = {"t": float(s["t"])}
            for n in ("phi", "E", "rho"):
                if n in s:
                    d[n] = np.asarray(s[n])
            for sp in ("e", "i"):
                if f"f_{sp}" in s:
                    m = s.meta.get(f"mass_{sp}", 1.0 if sp == "e" else 1836.0)
                    for q, a in zip("nuTq", moments(np.asarray(s[f"f_{sp}"], float), np.asarray(s[f"v_{sp}"]), m,
                                                    dv=dv_of(s, sp))):
                        d[f"{q}_{sp}"] = a
        return d

    def var(self, d, name):
        if name in d:
            return d[name]
        if name == "dn":
            return d["n_i"] - d["n_e"]
        if name == "lamD":
            return np.sqrt(d["T_e"] / np.maximum(d["n_e"], 1e-30))
        if name == "dxlam":
            return self.dx / self.lam / np.sqrt(d["T_e"] / np.maximum(d["n_e"], 1e-30))
        raise KeyError(f"{self.name}: 没有变量 {name}")

    # ---- 按壁面对折 ----
    def fold(self, a, name, side):
        """返回 (离壁距离 λ, 值)；side=full 时返回 (x, 值)。"""
        if side == "full":
            return self.x, a
        n = len(self.x) // 2
        d = (1.0 - np.abs(self.x[-n:])) / self.lam  # 右半 (由壁面向内递减) → 反序后由壁面向内
        d = d[::-1]
        left, right = a[:n], a[-n:][::-1]
        if name in VECTOR:
            left = -left  # 左壁指向壁面为 −x
        val = {"left": left, "right": right, "avg": 0.5 * (left + right)}[side]
        return d, val

    def fold_edges(self, side):
        """对折后各格的边界 (离壁距离 λ)。"""
        xf =np.concatenate([[self.x[0] - 0.5 * self.dx[0]], self.x + 0.5 * self.dx])
        if side == "full":
            return xf
        n = len(self.x) // 2
        return (1.0 - xf[-n - 1:])[::-1] / self.lam

    def interfaces(self):
        """格宽跳变处 (离壁距离 λ)：refine / stretch 网格的粗细界面。"""
        n = len(self.x) // 2
        h = self.dx[-n:][::-1]
        e = np.concatenate([[0.0], np.cumsum(h)]) / self.lam
        jump = np.where(np.abs(h[1:] / h[:-1] - 1) > 0.5)[0]
        return e[jump + 1]


def tlist(s):
    return [float(v) for v in s.split(",")] if s else []


def save(fig, args, default):
    out = args.out or default
    fig.savefig(out, dpi=args.dpi, bbox_inches="tight")
    print("图已保存:", out)
    if args.show:
        plt.show()


def setup_daxis(ax, args, src):
    if args.side == "full":
        ax.set_xlabel("x [L]")
        return
    ax.set_xlabel("distance from wall [λ]")
    if args.logd:
        ax.set_xscale("log")
    if args.dmax:
        ax.set_xlim(right=args.dmax)
    for d in src.interfaces():
        if not args.dmax or d < args.dmax:
            ax.axvline(d, color="#8a8984", lw=0.8, ls=(0, (2, 2)), zorder=0)


# ======================= 子命令 =======================
def cmd_prof(args):
    srcs = [Source(r) for r in args.runs]
    names = args.vars.split(",")
    times = tlist(args.t) or list(np.linspace(srcs[0].t[0], srcs[0].t[-1], 5))
    nc = min(len(names), args.cols)
    nr = int(np.ceil(len(names) / nc))
    fig, axs = plt.subplots(nr, nc, figsize=(4.2 * nc, 2.9 * nr), squeeze=False)
    multi = len(srcs) > 1
    for si, src in enumerate(srcs):
        for ti, t in enumerate(times):
            k = src.nearest(t)
            d = src.raw(k)
            if multi:  # 颜色 = 算例，线型 = 时刻
                color, ls = CAT[si % len(CAT)], DASH[ti % len(DASH)]
            else:  # 颜色 = 时刻 (蓝阶，由浅到深)
                color, ls = BLUE[int(round(ti * (len(BLUE) - 1) / max(len(times) - 1, 1)))], "-"
            lab = (f"{src.name} " if multi else "") + f"t={src.t[k]:.3g}"
            for vi, name in enumerate(names):
                ax = axs[vi // nc][vi % nc]
                try:
                    a = src.var(d, name)
                except KeyError as e:
                    print(e)
                    continue
                xx, yy = src.fold(np.asarray(a), name, args.side)
                ax.plot(xx, yy, color=color, ls=ls, label=lab, drawstyle="steps-mid" if args.steps else "default")
    for vi, name in enumerate(names):
        ax = axs[vi // nc][vi % nc]
        ax.set_ylabel(LABEL.get(name, name))
        setup_daxis(ax, args, srcs[0])
        if name in SIGNED:
            ax.axhline(0, color="#8a8984", lw=0.6)
        if args.logy and name not in SIGNED:
            ax.set_yscale("log")
    for vi in range(len(names), nr * nc):
        axs[vi // nc][vi % nc].set_visible(False)
    seen = {}
    for ax in axs.flat:
        for h, lab in zip(*ax.get_legend_handles_labels()):
            seen.setdefault(lab, h)
    axs[0][0].legend(seen.values(), seen.keys(), fontsize=7, ncol=1 + (len(seen) > 8))
    fig.suptitle(" vs ".join(s.name for s in srcs) + f"  ({args.side})", fontsize=10)
    fig.tight_layout()
    save(fig, args, os.path.join(srcs[0].path, f"prof_{'_'.join(names)}.png"))


def cmd_xt(args):
    src = Source(args.runs[0])
    ref = Source(args.ref) if args.ref else None
    ks = np.arange(len(src.t))
    if args.tmin is not None:
        ks = ks[src.t[ks] >= args.tmin]
    if args.tmax is not None:
        ks = ks[src.t[ks] <= args.tmax]
    rows = []
    for k in ks:
        a = np.asarray(src.var(src.raw(k), args.var))
        if ref is not None:
            kr = ref.nearest(src.t[k])
            a = a - np.asarray(ref.var(ref.raw(kr), args.var))
        rows.append(src.fold(a, args.var, args.side)[1])
    Z = np.array(rows)
    t = src.t[ks]
    te = np.concatenate([[t[0] - 0.5 * (t[1] - t[0])], 0.5 * (t[1:] + t[:-1]), [t[-1] + 0.5 * (t[-1] - t[-2])]])
    de = src.fold_edges(args.side)
    if args.logd and args.side != "full":
        de = de.copy()
        de[0] = 0.5 * de[1]
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    signed = args.var in SIGNED or ref is not None
    if signed:
        lim = args.vmax or np.nanpercentile(np.abs(Z), 99.5) or 1.0
        norm = TwoSlopeNorm(0.0, -lim, lim)
        cmap = DIV
    elif args.logz:
        pos = Z[Z > 0]
        norm, cmap = LogNorm(vmin=args.vmin or pos.min(), vmax=args.vmax or pos.max()), SEQ
    else:
        norm, cmap = plt.Normalize(args.vmin if args.vmin is not None else Z.min(), args.vmax or Z.max()), SEQ
    pc = ax.pcolormesh(de, te, Z, norm=norm, cmap=cmap, shading="flat", rasterized=True)
    cb = fig.colorbar(pc, ax=ax)
    cb.set_label(LABEL.get(args.var, args.var) + (f" − {ref.name}" if ref else ""))
    ax.set_ylabel("t [τi]")
    ax.grid(False)
    setup_daxis(ax, args, src)
    ax.set_title(f"{src.name}: {args.var}" + (f" − {ref.name}" if ref else "") + f" ({args.side})")
    fig.tight_layout()
    save(fig, args, os.path.join(src.path, f"xt_{args.var}{'_rel' if ref else ''}.png"))


def cmd_f(args):
    run = Run(args.runs[0])
    if not run.snap_files:
        sys.exit("没有分布函数快照")
    s = run.snap(t=float(args.t or 0))
    f = np.asarray(s[f"f_{args.sp}"], float)
    v = np.asarray(s[f"v_{args.sp}"])
    lam = s.meta.get("lambda", run.case.get("lambda", 1e-3))
    x = np.asarray(s["x"])
    xf = np.asarray(s["x_face"]) if "x_face" in s else None
    if xf is None:
        xf = np.concatenate([[x[0] - 0.5 * (x[1] - x[0])], 0.5 * (x[1:] + x[:-1]), [x[-1] + 0.5 * (x[-1] - x[-2])]])
    ve = v_edges(s, args.sp)  # 非均匀 v 网格时为真实格面 (速度网格关于 0 对称，翻转 v 时格面不变)
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    if args.side == "full":
        X, F = xf, f
        ax.set_xlabel("x [L]")
    else:  # 只画一侧壁面 (left / right，avg 视为 right)；v 取指向壁面为正
        n = len(x) // 2
        if args.side == "left":
            X, F = (xf[: n + 1] + 1.0) / lam, f[:n][:, ::-1]
        else:
            X, F = (1.0 - xf[-n - 1:])[::-1] / lam, f[-n:][::-1]
        ax.set_xlabel("distance from wall [λ]")
        if args.logd:
            X = X.copy()
            X[0] = 0.5 * X[1]
            ax.set_xscale("log")
        if args.dmax:
            ax.set_xlim(right=args.dmax)
    Z = F.T
    if args.log:
        pos = Z[Z > 0]
        vmin = args.vmin or max(pos.min(), pos.max() * 1e-8)
        norm = LogNorm(vmin=vmin, vmax=pos.max())
        Z = np.where(Z > vmin, Z, vmin)
    else:
        norm = plt.Normalize(0, Z.max())
    pc = ax.pcolormesh(X, ve, Z, norm=norm, cmap=SEQ, shading="flat", rasterized=True)
    fig.colorbar(pc, ax=ax).set_label(f"f_{args.sp} [n0/v_Te]")
    ax.set_ylabel(f"v{'' if args.side == 'full' else ' (toward wall +)'} [v_Te]")
    vlim = args.vlim
    if not vlim:  # 自动: f > 1e-6·max 的速度范围再放宽 20%
        hit = np.abs(v)[(f > 1e-6 * f.max()).any(axis=0)]
        vlim = 1.2 * hit.max() if hit.size else None
    if vlim:
        ax.set_ylim(-vlim, vlim)
    ax.grid(False)
    ax.set_title(f"{os.path.basename(run.path.rstrip('/'))}: f_{args.sp}, t = {float(s['t']):.3g} τi ({args.side})")
    fig.tight_layout()
    save(fig, args, os.path.join(run.path, f"f_{args.sp}_t{float(s['t']):.3g}.png"))


def cmd_fv(args):
    ds = tlist(args.d) or [0.25, 10, 100]
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    runs = [Run(r) for r in args.runs]
    multi = len(runs) > 1
    for ri, run in enumerate(runs):
        s = run.snap(t=float(args.t or 0))
        f = np.asarray(s[f"f_{args.sp}"], float)
        v = np.asarray(s[f"v_{args.sp}"])
        x = np.asarray(s["x"])
        lam = s.meta.get("lambda", run.case.get("lambda", 1e-3))
        for di, d in enumerate(ds):
            if args.side == "left":
                i = int(np.argmin(np.abs((x + 1) / lam - d)))
                fv, sgn = f[i], -1
            else:
                i = int(np.argmin(np.abs((1 - x) / lam - d)))
                fv, sgn = f[i], 1
            vv = sgn * v
            o = np.argsort(vv)
            color = CAT[ri % len(CAT)] if multi else CAT[di % len(CAT)]
            ls = DASH[di % len(DASH)] if multi else "-"
            dd = (1 - abs(x[i])) / lam
            ax.plot(vv[o], fv[o], color=color, ls=ls,
                    label=(f"{os.path.basename(run.path.rstrip('/'))} " if multi else "") + f"d={dd:.3g}λ, t={float(s['t']):.3g}")
    ax.set_xlabel("v (toward wall +) [v_Te]")
    ax.set_ylabel(f"f_{args.sp} [n0/v_Te]")
    if args.log:
        ax.set_yscale("log")
    if args.vlim:
        ax.set_xlim(-args.vlim, args.vlim)
    ax.legend(fontsize=7)
    fig.tight_layout()
    save(fig, args, os.path.join(runs[0].path, f"fv_{args.sp}.png"))


def cmd_diag(args):
    runs = [Run(r) for r in args.runs]
    cols = args.cols.split(",")
    fig, axs = plt.subplots(len(cols), 1, figsize=(6.4, 2.4 * len(cols)), squeeze=False, sharex=True)
    d0 = runs[0].diag
    for ri, run in enumerate(runs):
        d = run.diag
        for ci, c in enumerate(cols):
            ax = axs[ci][0]
            if c not in d:
                print(f"{run.path}: 没有列 {c} (可用: {', '.join(k for k in d if k != 't')})")
                continue
            y = d[c]
            if args.rel and ri > 0:
                y = y / np.interp(d["t"], d0["t"], d0[c]) - 1
            elif args.rel:
                continue
            ax.plot(d["t"], y, color=CAT[ri % len(CAT)], label=os.path.basename(run.path.rstrip("/")))
    for ci, c in enumerate(cols):
        ax = axs[ci][0]
        ax.set_ylabel(c + (f" / {os.path.basename(runs[0].path.rstrip('/'))} − 1" if args.rel else ""))
        if args.logy and not args.rel:
            ax.set_yscale("log")
        if args.tmax:
            ax.set_xlim(right=args.tmax)
    axs[-1][0].set_xlabel("t [τi]")
    axs[0][0].legend(fontsize=7)
    fig.tight_layout()
    save(fig, args, os.path.join(runs[0].path, f"diag{'_rel' if args.rel else ''}.png"))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["prof", "xt", "f", "fv", "diag"])
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--t", help="时刻 (τi)，逗号分隔")
    ap.add_argument("--vars", default="E,phi,dn,n_e,T_e,dxlam", help="prof 的变量")
    ap.add_argument("--var", default="E", help="xt 的变量")
    ap.add_argument("--ref", help="xt: 减去该对照算例 (同一网格)")
    ap.add_argument("--cols", default="ElectricEnergy,Tw.e-,flux.H+", help="diag 的列；prof 时为每行子图数")
    ap.add_argument("--sp", default="e", choices=["e", "i"])
    ap.add_argument("--d", help="fv: 离壁距离 (λ)，逗号分隔")
    ap.add_argument("--side", default="avg", choices=["avg", "left", "right", "full"])
    ap.add_argument("--dmax", type=float)
    ap.add_argument("--logd", action="store_true")
    ap.add_argument("--logy", action="store_true")
    ap.add_argument("--logz", action="store_true")
    ap.add_argument("--log", action="store_true", help="f / fv 用对数色标或对数轴")
    ap.add_argument("--steps", action="store_true", help="prof 画成格平均台阶")
    ap.add_argument("--rel", action="store_true")
    ap.add_argument("--tmin", type=float)
    ap.add_argument("--tmax", type=float)
    ap.add_argument("--vmin", type=float)
    ap.add_argument("--vmax", type=float)
    ap.add_argument("--vlim", type=float, help="速度轴范围 ±vlim")
    ap.add_argument("--out")
    ap.add_argument("--dpi", type=int, default=150)
    ap.add_argument("--show", action="store_true")
    args = ap.parse_args()
    if args.cmd == "prof":
        args.cols = int(args.cols) if args.cols.isdigit() else 3
    {"prof": cmd_prof, "xt": cmd_xt, "f": cmd_f, "fv": cmd_fv, "diag": cmd_diag}[args.cmd](args)


if __name__ == "__main__":
    main()
