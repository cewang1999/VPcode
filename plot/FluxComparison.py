"""
靶板粒子流 / 能流：离子、电子与离子自由流解析解对比 (Manfredi 2011, Fig.1)
用法: python3 plot/FluxComparison.py [runtime_diags.csv] [运行目录(含 results/ Grid/ case/, 可选)]
依赖: numpy, matplotlib
单位: 时间 tau_i, 粒子流 n_ELM*v_Ti, 能流 n_ELM*v_Ti*T_ELM
"""
import sys
import csv
import glob
import os

try:
    import numpy as np
    import matplotlib.pyplot as plt
except ImportError as err:
    sys.exit(f"缺少依赖: {err.name}。请用装有 numpy 和 matplotlib 的环境运行，例如:\n"
             "  ~/plotenv/bin/python plot/FluxComparison.py")

M_H = 1836.0
SIGMA = 0.1  # 源的空间宽度 (L_par)
OUT_EVERY = 100  # main.cpp 中输出 f_e 的步数间隔

ION, ELEC, REF = "#2a78d6", "#eb6834", "#8a8984"
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"


def read_diags(path):
    if not os.path.exists(path):
        sys.exit(f"找不到 {path}：请先在该目录运行程序生成 runtime_diags.csv，或在命令行指定其路径")
    with open(path) as fin:
        rows = list(csv.reader(fin))
    head, body = rows[0], np.array(rows[1:], dtype=float)
    return {name: body[:, k] for k, name in enumerate(head)}


def read_dt(run_dir, case="2"):
    """从 case 文件读取 dt (电子单位)"""
    fn = os.path.join(run_dir, "case", f"Case{case}.dat")
    if os.path.exists(fn):
        for line in open(fn):
            item = line.split()
            if len(item) >= 2 and item[0] == "dt":
                return float(item[1])
    return 0.02


def free_streaming(t, x=1.0, sigma=SIGMA):
    """高斯空间源 + 瞬时源的离子自由流解 (离子单位, m=1, T=1)"""
    s2 = sigma**2 + t**2
    w2 = sigma**2 / s2  # 速度方差
    u = x * t / s2  # 平均速度
    n = np.sqrt(w2) * np.exp(-0.5 * x**2 / s2)
    gamma = n * u
    q = n * (0.5 * (u**3 + 3 * u * w2) + u)  # Q_par + Q_perp(=Gamma*T)
    return gamma, q


def one_way_electron_flux(run_dir, dt):
    """旧诊断: 左壁第一个网格中心 v<0 部分的单向电子流"""
    grid = os.path.join(run_dir, "Grid", "Grid.txt")
    files = sorted(glob.glob(os.path.join(run_dir, "results", "vlasov_*.txt")))
    if not os.path.exists(grid) or len(files) < 2:
        return None
    v = np.loadtxt(grid, delimiter=",", skiprows=1)
    dv = v[1] - v[0]
    t, g = [], []
    for fn in files:
        step = int(os.path.basename(fn)[7:13])
        if step == 0:
            continue
        f0 = np.loadtxt(fn, delimiter=",", usecols=0)  # 第 0 列 = 左壁网格
        t.append(step * dt / np.sqrt(M_H))
        g.append(-np.sum(f0[v < 0] * v[v < 0]) * dv * np.sqrt(M_H))
    return np.array(t), np.array(g)


def style(ax, ylabel, tmax):
    ax.set_xlabel(r"time ($\tau_{\parallel i}$)", color=MUTED)
    ax.set_ylabel(ylabel, color=MUTED)
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=MUTED)
    ax.set_xlim(0, tmax)


def main():
    csv_file = sys.argv[1] if len(sys.argv) > 1 else "runtime_diags.csv"
    run_dir = sys.argv[2] if len(sys.argv) > 2 else os.path.dirname(os.path.abspath(csv_file))
    d = read_diags(csv_file)
    t = d["ts"]
    tmax = max(t[-1], 0.5)
    tf = np.linspace(1e-3, tmax, 800)
    gf, qf = free_streaming(tf)

    fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(15, 4.2))

    a1.plot(tf, gf, ":", color=REF, lw=2, label="ion, free-streaming")
    a1.plot(t, d["flux.H+"], "-", color=ION, lw=2, label="ion")
    a1.plot(t, d["flux.e-"], "--", color=ELEC, lw=2, label="electron (net, wall face)")
    old = one_way_electron_flux(run_dir, read_dt(run_dir))
    if old is not None:
        a1.plot(*old, "o", mfc="none", mec=ELEC, ms=8, mew=1.5, label="electron, old diagnostic (one-way)")
    style(a1, r"$\Gamma$ ($n_{ELM} v_{Ti}$)", tmax)
    a1.set_ylim(bottom=0)
    a1.set_title("Particle flux", color=INK, loc="left")
    a1.legend(frameon=False, labelcolor=INK)

    a2.plot(tf, qf, ":", color=REF, lw=2, label="ion, free-streaming")
    a2.plot(t, d["Enflux.H+"], "-", color=ION, lw=2, label="ion")
    a2.plot(t, d["Enflux.e-"], "--", color=ELEC, lw=2, label="electron")
    style(a2, r"$Q$ ($n_{ELM} v_{Ti} T_{ELM}$)", tmax)
    a2.set_ylim(bottom=0)
    a2.set_title("Energy flux", color=INK, loc="left")
    a2.legend(frameon=False, labelcolor=INK)

    # 电子流相对离子流的偏差，只画离子流不太小的时段
    ok = d["flux.H+"] > 0.05 * d["flux.H+"].max()
    a3.axhline(0, color=MUTED, lw=1)
    a3.plot(t[ok], 100 * (d["flux.e-"][ok] / d["flux.H+"][ok] - 1), "-", color=ELEC, lw=2)
    style(a3, r"$\Gamma_e/\Gamma_i - 1$ (%)", tmax)
    a3.set_title(r"Electron excess over ions", color=INK, loc="left")

    fig.tight_layout()
    out = os.path.join(os.path.dirname(os.path.abspath(csv_file)), "FluxComparison.png")
    fig.savefig(out, dpi=160)
    print("saved", out)
    plt.show()


if __name__ == "__main__":
    main()
