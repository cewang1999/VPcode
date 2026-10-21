"""
靶板粒子流 / 能流：离子、电子与离子自由流解析解对比 (Manfredi 2011, Fig.1)
用法: python FluxComparison.py [runtime_diags.csv] [运行目录(含 results/ 与 Grid/, 可选)]
单位: 时间 tau_i, 粒子流 n_ELM*v_Ti, 能流 n_ELM*v_Ti*T_ELM
"""
import sys
import glob
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

M_H = 1836.0
SIGMA = 0.1  # 源的空间宽度 (L_par)
DT = 0.02  # 电子单位时间步, 与 case 文件一致
OUT_EVERY = 100  # main.cpp 中输出 f_e 的步数间隔

ION, ELEC, REF = "#2a78d6", "#eb6834", "#8a8984"
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"


def free_streaming(t, x=1.0, sigma=SIGMA):
    """高斯空间源 + 瞬时源的离子自由流解 (离子单位, m=1, T=1)"""
    s2 = sigma**2 + t**2
    w2 = sigma**2 / s2  # 速度方差
    u = x * t / s2  # 平均速度
    n = np.sqrt(w2) * np.exp(-0.5 * x**2 / s2)
    gamma = n * u
    q = n * (0.5 * (u**3 + 3 * u * w2) + u)  # Q_par + Q_perp(=Gamma*T)
    return gamma, q


def one_way_electron_flux(run_dir):
    """旧诊断: 左壁第一个网格中心 v<0 部分的单向电子流"""
    grid = os.path.join(run_dir, "Grid", "Grid.txt")
    files = sorted(glob.glob(os.path.join(run_dir, "results", "vlasov_*.txt")))
    if not os.path.exists(grid) or not files:
        return None
    v = np.loadtxt(grid, delimiter=",", skiprows=1)
    dv = v[1] - v[0]
    t, g = [], []
    for fn in files:
        step = int(os.path.basename(fn)[7:13])
        if step == 0:
            continue
        f0 = np.loadtxt(fn, delimiter=",", usecols=0)  # 第 0 列 = 左壁网格
        t.append(step * DT / np.sqrt(M_H))
        g.append(-np.sum(f0[v < 0] * v[v < 0]) * dv * np.sqrt(M_H))
    return np.array(t), np.array(g)


def style(ax, ylabel):
    ax.set_xlabel(r"time ($\tau_{\parallel i}$)", color=MUTED)
    ax.set_ylabel(ylabel, color=MUTED)
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=MUTED)
    ax.set_xlim(0, 4)
    ax.set_ylim(bottom=0)


def main():
    csv = sys.argv[1] if len(sys.argv) > 1 else "runtime_diags.csv"
    run_dir = sys.argv[2] if len(sys.argv) > 2 else os.path.dirname(os.path.abspath(csv))
    d = pd.read_csv(csv)
    t = d["ts"].values
    tf = np.linspace(1e-3, 4, 800)
    gf, qf = free_streaming(tf)

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2))

    a1.plot(tf, gf, ":", color=REF, lw=2, label="ion, free-streaming")
    a1.plot(t, d["flux.H+"], "-", color=ION, lw=2, label="ion")
    a1.plot(t, d["flux.e-"], "--", color=ELEC, lw=2, label="electron (net, wall face)")
    old = one_way_electron_flux(run_dir)
    if old is not None:
        a1.plot(*old, "o", mfc="none", mec=ELEC, ms=8, mew=1.5, label="electron, old diagnostic (one-way)")
    style(a1, r"$\Gamma$ ($n_{ELM} v_{Ti}$)")
    a1.set_title("Particle flux", color=INK, loc="left")
    a1.legend(frameon=False, labelcolor=INK)

    a2.plot(tf, qf, ":", color=REF, lw=2, label="ion, free-streaming")
    a2.plot(t, d["Enflux.H+"], "-", color=ION, lw=2, label="ion")
    a2.plot(t, d["Enflux.e-"], "--", color=ELEC, lw=2, label="electron")
    style(a2, r"$Q$ ($n_{ELM} v_{Ti} T_{ELM}$)")
    a2.set_title("Energy flux", color=INK, loc="left")
    a2.legend(frameon=False, labelcolor=INK)

    fig.tight_layout()
    out = os.path.join(os.path.dirname(os.path.abspath(csv)), "FluxComparison.png")
    fig.savefig(out, dpi=160)
    print("saved", out)


if __name__ == "__main__":
    main()
