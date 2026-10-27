"""
kobra 诊断 (quasi-kobra-sup 中 usr.c 末尾的 diag_mon 输出 diag_budget.dat) 与 VPcode 逐项对比
  通量: 净通量 Γ = -dN/dt / 2 (两壁平均) 与贴壁格出射值 × v (kobra 的 wallinfo / VPcode 的 flux1)
  Γe/Γi - 1, 早期电势能峰 W/E0×1000, 高斯定律 Q = E_R - E_L, 壁面 / 全域电子温度, 离子 T_i·t²
用法: ~/plotenv/bin/python plot/KobraDiag.py kobra数据目录 [VPcode运行目录 ...]
      kobra 数据目录可用 tools/fetch_kobra.sh 拷到 runs/kobra_<标签>/
单位换算: kobra 中 λ = 1、密度乘 S (由 N_i(0) 自动求得, 通常 1e6)、时间为 L/v_Te；
          这里统一换成 VPcode 的 n_elm、n v_Ti、τi
"""
import os
import sys

try:
    import numpy as np
    import matplotlib.pyplot as plt
except ImportError as err:
    sys.exit(f"缺少依赖: {err.name}。请用 ~/plotenv/bin/python 运行")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import PaperFigs as PF  # noqa: E402  (read_diags, run_label, paper_axes, 字体设置)

N0_UNIT = np.sqrt(2 * np.pi) * PF.SIGMA  # 初始高斯 (σ = 0.1) 的总粒子数 (n_elm L)


def read_kobra(kdir):
    path = os.path.join(kdir, "diag_budget.dat")
    if not os.path.exists(path):
        sys.exit(f"找不到 {path}")
    head = [ln for ln in open(path) if ln.startswith("#")]
    names = head[2][1:].split()
    par = head[1][1:].split()
    mass = [float(x) for x in par[par.index("m") + 1:par.index("alpha")]]
    body = np.loadtxt(path)
    if body.ndim == 1:
        body = body[None, :]
    c = {k: body[:, i] for i, k in enumerate(names)}
    mi = mass[1] / mass[0]
    S = c["N1"][0] / N0_UNIT
    k = {"S": S, "mi": mi, "ts": c["t"] / np.sqrt(mi)}
    t = c["t"]
    for s, sp in ((0, "e-"), (1, "H+")):
        N = c[f"N{s}"]
        k["N." + sp] = N / S
        k["out." + sp] = (N[0] - N) / S / 2  # 单壁累计流出 (n L)
        k["flux." + sp] = -np.gradient(N, t) / 2 / S * np.sqrt(mi)  # n v_Ti
        k["flux1." + sp] = 0.5 * (c[f"L_j{s}"] + c[f"R_j{s}"]) / S * np.sqrt(mi)
        k["K." + sp] = c[f"K{s}"] / S
        k["nw." + sp] = 0.5 * (c[f"L_n{s}"] + c[f"R_n{s}"]) / S
        k["Tw." + sp] = 0.5 * (c[f"L_T{s}"] + c[f"R_T{s}"])
        k["floor." + sp] = c[f"floormass{s}"] / S
    k["W"] = c["W"] / S  # λ²/2 ΣE²dx (n T L)
    k["Q"] = c["Q"] / S
    k["dE"] = (c["ER"] - c["EL"]) / S  # λ²(E_R - E_L)
    k["gmax"] = c["gmax"] / np.maximum(c["rhomax"], 1e-30)
    return k


def from_vpcode(run_dir):
    d = PF.read_diags(run_dir)
    lam = 1e-3
    case_dir = os.path.join(run_dir, "case")
    if os.path.isdir(case_dir):
        for fn in os.listdir(case_dir):
            for line in open(os.path.join(case_dir, fn)):
                item = line.split()
                if len(item) >= 2 and item[0] == "lambda":
                    lam = float(item[1])
    t = d["ts"]
    v = {"ts": t, "W": 0.5 * lam**2 * d["ElectricEnergy"]}
    for sp in ("e-", "H+"):
        v["flux." + sp] = d["flux." + sp]
        v["flux1." + sp] = d["flux1." + sp]
        v["N." + sp] = d["N." + sp]
        v["K." + sp] = d["totalKE." + sp]
        v["nw." + sp] = d["nw." + sp]
        v["Tw." + sp] = d["Tw." + sp]
        g = d["flux." + sp]
        v["out." + sp] = np.concatenate([[0], np.cumsum(0.5 * (g[1:] + g[:-1]) * np.diff(t))])
    return v


def cum(y, t):
    return np.concatenate([[0], np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(t))])


def e0(r):
    return r["K.e-"][0] + r["K.H+"][0] + r["N.e-"][0] + r["N.H+"][0]  # 平行动能 + 垂直能量 N·T


def summary(name, r, kobra):
    t = r["ts"]
    W = r["W"] / e0(r) * 1000
    k = np.argmax(np.where(t < 0.05, W, 0))
    print(f"== {name}: 早期电势能峰 W/E0×1000 = {W[k]:.3f} (t = {t[k]:.4f} τi)")
    cell_e, cell_i = cum(r["flux1.e-"], t), cum(r["flux1.H+"], t)
    for x in (0.5, 1, 2, 3, 4):
        if t[-1] < x - 0.01:
            continue
        j = np.argmin(abs(t - x))
        net = 100 * (r["out.e-"][j] / r["out.H+"][j] - 1)
        cel = 100 * (cell_e[j] / cell_i[j] - 1)
        Te = 2 * r["K.e-"][j] / r["N.e-"][j]
        line = (f"  t={t[j]:4.2f}: 累积Γe/Γi-1 净 {net:+6.2f}%  贴壁格 {cel:+6.2f}%  "
                f"壁面Te {r['Tw.e-'][j]:.2e}  全域Te {Te:.2e}  Ti·t² {r['Tw.H+'][j] * t[j]**2:.4f}")
        if kobra:
            err = (r["Q"][j] - r["dE"][j]) / max(abs(r["Q"]).max(), 1e-30)
            line += f"  高斯误差 {100 * err:+.2f}%  floor {r['floor.e-'][j]:.1e}"
        print(line)


def plot(kdir, k, vps):
    fig, ax = plt.subplots(2, 3, figsize=(15, 8))
    ax = ax.ravel()
    runs = [("kobra", k, "k")] + [(os.path.basename(os.path.normpath(p)), v, c)
                                  for (p, v), c in zip(vps, ["tab:red", "tab:blue", "tab:green"])]
    for name, r, col in runs:
        t = r["ts"]
        ax[0].plot(t, r["flux.H+"], "-", color=col, lw=1.2, label=f"{name} ion net")
        ax[0].plot(t, r["flux.e-"], "--", color=col, lw=1.2, label=f"{name} e net")
        ax[0].plot(t, r["flux1.e-"], ":", color=col, lw=1.0, label=f"{name} e cell")
        ok = r["flux.H+"] > 0.05 * r["flux.H+"].max()
        ax[1].plot(t[ok], 100 * (r["flux.e-"][ok] / r["flux.H+"][ok] - 1), "-", color=col, lw=1.0,
                   label=f"{name} net")
        ax[1].plot(t[ok], 100 * (r["flux1.e-"][ok] / r["flux1.H+"][ok] - 1), ":", color=col, lw=1.0,
                   label=f"{name} cell")
        pos = t > 0
        ax[2].semilogx(t[pos], (r["W"] / e0(r) * 1000)[pos], "-", color=col, lw=1.2, label=name)
        ax[4].semilogy(t, np.maximum(r["Tw.e-"], 1e-6), "-", color=col, lw=1.2, label=f"{name} wall")
        ax[4].semilogy(t, 2 * r["K.e-"] / r["N.e-"], "--", color=col, lw=1.0, label=f"{name} global")
        late = t > 0.2
        ax[5].plot(t[late], (r["Tw.H+"] * t**2)[late], "-", color=col, lw=1.2, label=name)
    t = k["ts"]
    ax[3].plot(t, 100 * (k["Q"] - k["dE"]) / max(abs(k["Q"]).max(), 1e-30), "k-", lw=1.0,
               label=r"$(Q-\lambda^2(E_R-E_L))/\max|Q|$")
    ax[3].plot(t, 100 * k["gmax"], "k:", lw=1.0, label=r"max$|g|$/max$|\rho|$")
    ax[5].axhline(PF.SIGMA**2 / 2, color="0.5", lw=0.8)
    titles = [r"$\Gamma$ ($n v_{Ti}$)", r"$\Gamma_e/\Gamma_i-1$ (%)", r"$W/E_0\times1000$",
              "kobra: Gauss / Poisson residual (%)", r"$T_e$ ($T_{ELM}$)", r"$T_i t^2$ (wall)"]
    for a, ti in zip(ax, titles):
        a.set_title(ti, fontsize=11)
        PF.paper_axes(a)
        a.legend(frameon=False, fontsize=7)
    ax[0].set_xlim(0, 4)
    ax[1].set_xlim(0, 4)
    ax[1].set_ylim(-5, 30)
    ax[2].set_xlim(1e-3, 4)
    fig.text(0.5, 0.005, "wall n, T: kobra from the outgoing part of the wall cell; VPcode from the wall face flux",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    out = os.path.join(kdir, "KobraDiag.png")
    fig.savefig(out, dpi=150)
    return out


def main():
    if len(sys.argv) < 2:
        sys.exit("用法: ~/plotenv/bin/python plot/KobraDiag.py kobra数据目录 [VPcode运行目录 ...]")
    kdir = sys.argv[1]
    k = read_kobra(kdir)
    print(f"kobra: 密度缩放 S = {k['S']:.4g}, m_i/m_e = {k['mi']:.6g}, 数据到 t = {k['ts'][-1]:.3f} τi")
    summary("kobra", k, True)
    vps = [(p, from_vpcode(p)) for p in sys.argv[2:]]
    for p, v in vps:
        summary(os.path.basename(os.path.normpath(p)), v, False)
    print("saved", plot(kdir, k, vps))


if __name__ == "__main__":
    main()
