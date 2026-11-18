"""kobra 算例的壁面粒子流 (读 RUN/flux_mpi_nrel100.out，kobra 的 WALLINFO 输出，每步一行)。

用法:
  python plot/KobraFlux.py RUN... [--lam 0.01] [--log] [--tmax 2] [--out 文件名] [--show]
例:
  python plot/KobraFlux.py ~/Vlasov/kobra-amr/naga
  python plot/KobraFlux.py ~/Vlasov/kobra-amr/runs/KA4 ~/Vlasov/kobra-nu/runs/KV1U --log

文件列: t  Γ_e  Γ_i  格数  最细层格数 (Γ 为左壁出射粒子流，kobra 单位: 密度以 1/λ² 计、速度以 v_Te 计)。
文件以追加方式写入，同一目录跑过多次时时间会从头开始；只取最后一次 (最后一个时间回零之后) 的数据。
画图量: Γ/(n0 v_Te)，n0 = 1/λ² (λ 由 out/prof.vpd 的 meta lambda 给出，没有时用 --lam)；t 以 τi = √(m_i/m_e) L/v_Te 计。
上图: Γ_e、Γ_i (实线电子、虚线离子，每个算例一种颜色)；下图: Γ_e/Γ_i (对数轴；离子到壁之前很大，鞘层建立后趋于 1)。
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

from VPplot import CAT  # noqa: E402

MI = 1836.0


def read_lambda(run, default):
    path = os.path.join(run, "out", "prof.vpd")
    if os.path.exists(path):
        with open(path, "rb") as fh:
            for raw in fh:
                line = raw.decode("utf-8", "replace").strip()
                if line == "end":
                    break
                parts = line.split()
                if len(parts) >= 3 and parts[0] == "meta" and parts[1] == "lambda":
                    return float(parts[2])
    if default is None:
        sys.exit(f"{run}: 没有 out/prof.vpd 的 meta lambda，请用 --lam 给出 λ")
    return default


def read_flux(run):
    path = os.path.join(run, "flux_mpi_nrel100.out")
    if not os.path.exists(path):
        sys.exit(f"没有 {path} (算例需在 .dat 中加 wallinfo)")
    d = np.loadtxt(path, ndmin=2)
    t = d[:, 0]
    reset = np.where(np.diff(t) < 0)[0]
    if len(reset):
        d = d[reset[-1] + 1:]
    return d


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--lam", type=float, default=None, help="λ (以 L 计)，prof.vpd 中没有时使用")
    ap.add_argument("--log", action="store_true", help="粒子流用对数轴")
    ap.add_argument("--tmax", type=float, default=None, help="时间上限 (τi)")
    ap.add_argument("--out", default=None)
    ap.add_argument("--show", action="store_true")
    a = ap.parse_args()

    tau = np.sqrt(MI)
    fig, (ax, axr) = plt.subplots(2, 1, figsize=(8, 6.5), sharex=True, gridspec_kw={"height_ratios": [2, 1]})
    for k, run in enumerate(a.runs):
        run = os.path.expanduser(run)
        lam = read_lambda(run, a.lam)
        n0 = 1.0 / lam**2
        d = read_flux(run)
        t = d[:, 0] / tau
        ge = d[:, 1] / n0
        gi = d[:, 2] / n0
        if a.tmax is not None:
            m = t <= a.tmax
            t, ge, gi = t[m], ge[m], gi[m]
        c = CAT[k % len(CAT)]
        name = os.path.basename(os.path.normpath(run))
        ax.plot(t, ge, color=c, lw=1.6, label=f"{name}  Γ_e")
        ax.plot(t, gi, color=c, lw=1.6, ls="--", label=f"{name}  Γ_i")
        with np.errstate(divide="ignore", invalid="ignore"):
            r = np.where(gi > 1e-12 * max(gi.max(), 1e-300), ge / gi, np.nan)
        axr.plot(t, r, color=c, lw=1.4, label=name)
        print(f"{name}: λ = {lam:g}，t = 0 … {t[-1]:.4f} τi ({len(t)} 步)，"
              f"末时刻 Γ_e = {ge[-1]:.5g}，Γ_i = {gi[-1]:.5g}，Γ_e/Γ_i = {ge[-1] / gi[-1]:.4f}")
    if a.log:
        ax.set_yscale("log")
    ax.set_ylabel("wall particle flux Γ / (n0 v_Te)")
    ax.legend(frameon=False, fontsize=8, ncol=2 if len(a.runs) > 1 else 1)
    ax.grid(alpha=0.3)
    axr.axhline(1.0, color="#8a8984", lw=0.8)
    axr.set_ylabel("Γ_e / Γ_i")
    axr.set_xlabel("t [τi]")
    axr.grid(alpha=0.3)
    axr.set_yscale("log")
    fig.suptitle("wall particle fluxes (left wall)")
    fig.tight_layout()
    out = a.out or os.path.join(os.path.expanduser(a.runs[0]), "flux.png")
    fig.savefig(out, dpi=150)
    print(f"图已保存: {out}")
    if a.show:
        plt.show()


if __name__ == "__main__":
    main()
