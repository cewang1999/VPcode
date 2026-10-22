"""
按 Manfredi 2011 的坐标轴重画结果，便于和文献逐图对比
  Fig1.png : 粒子流 / 能流 (文献图 1)，含离子自由流解析解
  Fig7.png : 壁面密度 / 平行温度 (文献图 7)，需要 runtime_diags.csv 含 nw、Tw 列
  Ratio.png: 净通量 Γe/Γi - 1
用法: python3 plot/PaperFigs.py 运行目录1 [运行目录2 ...]
      每个运行目录需含 runtime_diags.csv (与 case/ 子目录, 用于读取 λ 等参数做图例)
依赖: numpy, matplotlib  (例如 ~/plotenv/bin/python)
"""
import sys
import os
import csv

try:
    import numpy as np
    import matplotlib.pyplot as plt
except ImportError as err:
    sys.exit(f"缺少依赖: {err.name}。请用装有 numpy 和 matplotlib 的环境运行，例如:\n"
             "  ~/plotenv/bin/python plot/PaperFigs.py runs/Case4")

SIGMA = 0.1  # 源的空间宽度 (L_par)
LINES = ["-", "--", ":", "-."]
plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "dejavuserif", "font.size": 12})


def read_diags(run_dir):
    path = os.path.join(run_dir, "runtime_diags.csv")
    if not os.path.exists(path):
        sys.exit(f"找不到 {path}")
    with open(path) as fin:
        rows = list(csv.reader(fin))
    head = rows[0]
    body = np.array([r for r in rows[1:] if len(r) == len(head)], dtype=float)
    return {name: body[:, k] for k, name in enumerate(head)}


def run_label(run_dir):
    """由 case 文件生成图例: Case 名, λ, Nx, Nv, dt"""
    name = os.path.basename(os.path.normpath(run_dir))
    par = {}
    case_dir = os.path.join(run_dir, "case")
    if os.path.isdir(case_dir):
        for fn in os.listdir(case_dir):
            for line in open(os.path.join(case_dir, fn)):
                item = line.split()
                if len(item) >= 2:
                    par[item[0]] = item[1]
    if not par:
        return name
    return (f"{name}: $\\lambda$={par.get('lambda', '1e-3')}, Nx={par.get('Nx')}, "
            f"Nv={par.get('Nv')}, dt={par.get('dt')}")


def free_streaming(t, x=1.0, sigma=SIGMA):
    """高斯空间源 + 瞬时源的离子自由流解 (离子单位): 密度, 粒子流, 能流, 平行温度"""
    s2 = sigma**2 + t**2
    w2 = sigma**2 / s2  # 速度方差 = 平行温度
    u = x * t / s2
    n = np.sqrt(w2) * np.exp(-0.5 * x**2 / s2)
    gamma = n * u
    q = n * (0.5 * (u**3 + 3 * u * w2) + u)  # Q_par + Q_perp(=Gamma*T)
    return n, gamma, q, w2


def paper_axes(ax, xlabel=r"time ($\tau_{//i}$)"):
    ax.set_xlabel(xlabel)
    ax.tick_params(direction="in", top=True, right=True, which="both")
    ax.minorticks_on()


def fig1(runs, out_dir):
    fig, axs = plt.subplots(len(runs), 2, figsize=(8.5, 3.8 * len(runs)), squeeze=False)
    tf = np.linspace(1e-3, 4, 800)
    _, gf, qf, _ = free_streaming(tf)
    for (run_dir, d), (a1, a2) in zip(runs, axs):
        t = d["ts"]
        a1.plot(t, d["flux.H+"], "k-", lw=1.4)
        a1.plot(t, d["flux.e-"], "k--", lw=1.4)
        a1.plot(tf, gf, "k:", lw=1.2)
        a2.plot(t, d["Enflux.H+"], "k-", lw=1.4)
        a2.plot(t, d["Enflux.e-"], "k--", lw=1.4)
        a2.plot(tf, qf, "k:", lw=1.2)
        a1.set_ylim(0, 0.15)
        a1.set_yticks([0, 0.05, 0.10, 0.15])
        a1.set_ylabel(r"$\Gamma$ ($n_{ELM}v_{Ti}$)")
        a2.set_ylim(0, 0.4)
        a2.set_yticks([0, 0.2, 0.4])
        a2.set_ylabel(r"Q ($n_{ELM}v_{Ti}T_{ELM}$)")
        for ax in (a1, a2):
            ax.set_xlim(0, 4)
            ax.set_xticks([0, 1, 2, 3, 4])
            paper_axes(ax)
        a1.set_title(run_label(run_dir), fontsize=9, loc="left")
    fig.text(0.5, 0.005, "solid: ions, dashed: electrons, dotted: ion free-streaming model",
             ha="center", fontsize=10)
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    fig.savefig(os.path.join(out_dir, "Fig1.png"), dpi=170)


def fig7(runs, out_dir):
    runs = [(r, d) for r, d in runs if "nw.e-" in d]
    if not runs:
        print("跳过 Fig7: runtime_diags.csv 中没有 nw / Tw 列（需用新版程序重新运行）")
        return False
    fig, axs = plt.subplots(len(runs), 2, figsize=(8.5, 3.8 * len(runs)), squeeze=False)
    tf = np.linspace(1e-3, 4, 800)
    nf, _, _, tif = free_streaming(tf)
    for (run_dir, d), (a1, a2) in zip(runs, axs):
        t = d["ts"]
        a1.plot(t, d["nw.H+"], "k-", lw=1.4)
        a1.plot(t, d["nw.e-"], "k--", lw=1.4)
        a1.plot(tf, nf, "k:", lw=1.2)
        a1.set_ylim(0, 0.08)
        a1.set_yticks([0, 0.02, 0.04, 0.06, 0.08])
        a1.set_ylabel(r"n / $n_{ELM}$")
        a1.set_title(run_label(run_dir), fontsize=9, loc="left")

        a2.plot(t, d["Tw.H+"], "k-", lw=1.4)
        a2.plot(t, d["Tw.e-"], "k--", lw=1.4)
        a2.plot(tf, tif, "k:", lw=1.0)
        a2.set_ylim(0, 0.010)
        a2.set_yticks(np.arange(0, 0.0101, 0.002))
        a2.set_ylabel(r"T / $T_{ELM}$")
        # 插图: 离子温度 × t^2 (文献插图)
        ins = a2.inset_axes([0.5, 0.5, 0.45, 0.42])
        ins.plot(t, d["Tw.H+"] * t**2, "k-", lw=1.0)
        ins.plot(tf, tif * tf**2, "k:", lw=0.8)
        ins.set_xlim(0, 4)
        ins.set_ylim(0, 0.003)
        ins.set_yticks([0, 0.001, 0.002, 0.003])
        ins.tick_params(labelsize=7, direction="in")
        ins.set_xlabel("time", fontsize=7)
        for ax in (a1, a2):
            ax.set_xlim(0, 4)
            ax.set_xticks([0, 1, 2, 3, 4])
            paper_axes(ax)
    fig.text(0.5, 0.005, "solid: ions, dashed: electrons, dotted: ion free-streaming model (eq. 11, Gaussian source)",
             ha="center", fontsize=10)
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    fig.savefig(os.path.join(out_dir, "Fig7.png"), dpi=170)
    return True


def ratio(runs, out_dir):
    fig, ax = plt.subplots(figsize=(7, 4))
    for k, (run_dir, d) in enumerate(runs):
        t = d["ts"]
        ok = d["flux.H+"] > 0.05 * d["flux.H+"].max()
        ax.plot(t[ok], 100 * (d["flux.e-"][ok] / d["flux.H+"][ok] - 1), "k" + LINES[k % len(LINES)],
                lw=1.2, label=run_label(run_dir))
    ax.axhline(0, color="0.5", lw=0.8)
    ax.set_xlim(0, 4)
    ax.set_ylabel(r"$\Gamma_e/\Gamma_i - 1$ (%)")
    paper_axes(ax)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "Ratio.png"), dpi=170)


def main():
    if len(sys.argv) < 2:
        sys.exit("用法: python3 plot/PaperFigs.py 运行目录1 [运行目录2 ...]")
    runs = [(r, read_diags(r)) for r in sys.argv[1:]]
    out_dir = runs[0][0] if len(runs) == 1 else os.path.dirname(os.path.abspath(os.path.normpath(runs[0][0])))
    fig1(runs, out_dir)
    got7 = fig7(runs, out_dir)
    ratio(runs, out_dir)
    print("saved Fig1.png" + (", Fig7.png" if got7 else "") + ", Ratio.png in", out_dir)
    plt.show()


if __name__ == "__main__":
    main()
