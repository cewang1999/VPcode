"""kobra (VPD 输出) 算例与 VPcode 算例的壁面密度 / 平行温度，按 Manfredi 2011 图 7 的坐标轴画在一起。

用法:
  python plot/KobraFig7.py RUN [RUN ...] [--tmax 4] [--out 文件名] [--show]
  有 runtime_diags.csv 的目录当作 VPcode 算例 (画线)，其余当作 kobra 算例 (画点)，顺序任意
例:
  python plot/KobraFig7.py runs/KB16 runs/KB16N256 runs/Case83 runs/Case87 runs/Case46

壁面量的定义同 VPcode 的 runtime_diags.csv (nw、Tw)：只取出射方向 (指向壁面) 的那一半分布，
  n_w = ∫f dv，T_w = m (<v²> − <v>²)，左右两壁平均。
  VPcode: 直接读 runtime_diags.csv (PFC3 用壁面面通量 / v，Scheme kobra 用贴壁格的 f)。
  kobra : 没有这个诊断，由快照 out/f_<ts>.vpd 的贴壁格 f 按同一定义计算 (每个快照一个点，画成圆点)。
上排: 文献的线性坐标 (n 0–0.08，T 0–0.010)；下排: 对数坐标 (后期 T_e 降到 1e-4 量级，线性轴上看不出差别)。
实线离子、虚线电子；点线为离子自由流解析解。
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import matplotlib  # noqa: E402

if "--show" not in sys.argv:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import PaperFigs as PF  # noqa: E402
from vpio import Run  # noqa: E402
from VPplot import CAT  # noqa: E402


def wall_moments(f, v, mass):
    """贴壁格的出射半分布的 n、T；f 为 (nx, nv)。左右两壁平均。"""
    dv = v[1] - v[0]
    ns, ts = [], []
    for row, m in ((f[0], v < 0), (f[-1], v > 0)):
        g, a = np.asarray(row, float)[m], np.abs(v[m])
        n = g.sum() * dv
        u = (g * a).sum() * dv / n
        ns.append(n)
        ts.append(mass * ((g * a * a).sum() * dv / n - u * u))
    return np.mean(ns), np.mean(ts)


def kobra_fixdt(path):
    """算例目录下 kobra .dat 中的 fixdt (外层步长，kobra 时间单位)；没有时返回 None。"""
    for fn in sorted(os.listdir(path)):
        if fn.endswith(".dat"):
            for ln in open(os.path.join(path, fn), errors="replace"):
                p = ln.split()
                if len(p) >= 2 and p[0].lower() == "fixdt":
                    return float(p[1])
    return None


def kobra_wall(path):
    """kobra 快照的壁面量。单精度 kobra 的时钟 (float 累加 dt) 有漂移 (2 τi 时 +0.15%，4 τi 时 +0.4%)，
    有 fixdt 时用 ts × fixdt 作为时间。"""
    run = Run(path)
    fixdt = kobra_fixdt(path)
    out = {k: [] for k in ("t", "nw.e-", "Tw.e-", "nw.H+", "Tw.H+")}
    for k in range(len(run.snap_files)):
        s = run.snap(k)
        nx = len(s["x"])
        if fixdt:
            out["t"].append(int(s["ts"]) * fixdt / np.sqrt(float(s.meta.get("mass_i", 1836.0))))
        else:
            out["t"].append(float(s["t"]))
        for sp, tag, mass in (("e", "e-", 1.0), ("i", "H+", 1836.0)):
            v = np.asarray(s[f"v_{sp}"])
            n, T = wall_moments(np.asarray(s[f"f_{sp}"]).reshape(nx, len(v)), v, float(s.meta.get(f"mass_{sp}", mass)))
            out[f"nw.{tag}"].append(n)
            out[f"Tw.{tag}"].append(T)
    return {k: np.array(v) for k, v in out.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+", help="算例目录；有 runtime_diags.csv 的当作 VPcode 算例 (画线)，否则当作 kobra 算例 (由快照算，画点)")
    ap.add_argument("--tmax", type=float, default=4.0, help="横轴上限 (τi)")
    ap.add_argument("--out", default=None)
    ap.add_argument("--show", action="store_true")
    a = ap.parse_args()

    runs = []  # (名字, 是否 kobra, 数据)
    for r in a.runs:
        r = os.path.expanduser(r)
        name = os.path.basename(os.path.normpath(r))
        if os.path.exists(os.path.join(r, "runtime_diags.csv")):
            d = PF.read_diags(r)
            d["t"] = d["ts"]
            runs.append((name, False, d))
        else:
            runs.append((name, True, kobra_wall(r)))
    if not any(k for _, k, _ in runs):
        print("注意: 没有 kobra 算例 (所有目录都有 runtime_diags.csv)")
    tf = np.linspace(1e-3, 4, 800)
    nf, _, _, tif = PF.free_streaming(tf)

    fig, axs = plt.subplots(2, 2, figsize=(10.5, 8.2))
    marks = ["o", "s", "^", "D", "v"]
    for row in range(2):
        a1, a2 = axs[row]
        nk = 0
        for k, (name, isk, d) in enumerate(runs):
            c = CAT[k % len(CAT)]
            if isk:
                mk = marks[nk % len(marks)]
                nk += 1
                a1.plot(d["t"], d["nw.H+"], mk, color=c, ms=4.5, label=name + " (snapshots)")
                a1.plot(d["t"], d["nw.e-"], mk, color=c, ms=4.5, mfc="none")
                a2.plot(d["t"], d["Tw.H+"], mk, color=c, ms=4.5, label=name + " (snapshots)")
                a2.plot(d["t"], d["Tw.e-"], mk, color=c, ms=4.5, mfc="none")
            else:
                a1.plot(d["t"], d["nw.H+"], color=c, lw=1.2, label=name)
                a1.plot(d["t"], d["nw.e-"], color=c, lw=1.2, ls="--")
                a2.plot(d["t"], d["Tw.H+"], color=c, lw=1.2, label=name)
                a2.plot(d["t"], d["Tw.e-"], color=c, lw=1.2, ls="--")
        a1.plot(tf, nf, "k:", lw=1.0)
        a2.plot(tf, tif, "k:", lw=1.0)
        a1.set_ylabel(r"n / $n_{ELM}$")
        a2.set_ylabel(r"T / $T_{ELM}$")
        if row == 0:
            a1.set_ylim(0, 0.08)
            a2.set_ylim(0, 0.010)
            a2.set_yticks(np.arange(0, 0.0101, 0.002))
            a1.set_title("wall density (axes of Manfredi 2011, Fig. 7)", fontsize=10)
            a2.set_title("wall parallel temperature", fontsize=10)
            a1.legend(frameon=False, fontsize=8)
        else:
            a1.set_yscale("log")
            a1.set_ylim(1e-3, 0.1)
            a2.set_yscale("log")
            a2.set_ylim(5e-5, 1.0)
        for ax in (a1, a2):
            ax.set_xlim(0, a.tmax)
            PF.paper_axes(ax)
    fig.text(0.5, 0.005, "solid / filled: ions, dashed / open: electrons, dotted: ion free-streaming model", ha="center", fontsize=10)
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    first_k = next((os.path.expanduser(r) for r, (_, isk, _) in zip(a.runs, runs) if isk), os.path.expanduser(a.runs[0]))
    out = a.out or os.path.join(first_k, "Fig7_cmp.png")
    fig.savefig(out, dpi=160)
    print(f"图已保存: {out}")

    # 数值表：各 kobra 算例快照时刻的壁面 T_e、T_i，及与各 VPcode 算例之比
    vp = [(n, d) for n, isk, d in runs if not isk]
    for name, isk, kw in runs:
        if not isk:
            continue
        print(f"\n{'t':>5} {name + ' Te_w':>18} {'Ti_w':>10} {'ne_w':>9} {'ni_w':>9} | " + " | ".join(f"Te/{n} Ti/.." for n, _ in vp))
        for i, t in enumerate(kw["t"]):
            cols = []
            for _, d in vp:
                if t > d["t"][-1] + 0.01:
                    cols.append(f"{'—':>16}")
                    continue
                j = int(np.argmin(np.abs(d["t"] - t)))
                cols.append(f"{kw['Tw.e-'][i] / d['Tw.e-'][j]:8.3f} {kw['Tw.H+'][i] / d['Tw.H+'][j]:7.3f}")
            print(f"{t:5.2f} {kw['Tw.e-'][i]:18.4e} {kw['Tw.H+'][i]:10.4e} {kw['nw.e-'][i]:9.5f} {kw['nw.H+'][i]:9.5f} | " + " | ".join(cols))
    if a.show:
        plt.show()


if __name__ == "__main__":
    main()
