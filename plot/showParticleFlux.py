import Freestreamingmodel as FS
import KobraOutput as KOP
import VPcode as VP
import numpy as np
import matplotlib.pyplot as plt


fsm = FS.FSM(
    FS.Setting(True, False, True, True), "Human-defined Model", FS.Ion("H+", 18.36, 1)
)
kobra = KOP.Kobra("/Users/wangce/Vlasov/VPcode/plot/flux.out")
vpcode = VP.VPcode("/Users/wangce/Vlasov/VPcode/runtime_diags.csv", 9)

# VPcode
t_vp = vpcode.getTime()
f_vp = vpcode.getCol() / np.sqrt(fsm.ions.mass)

# FSM
time = np.linspace(0.01, 40, 1000)
f_fsm = fsm.Output(time)
t_fsm = fsm.TimeTransfer(time)

# Kobra
t_k = kobra.get_column(0)
fe_k = kobra.get_column(1) / 1000000 / 2
fi_k = kobra.get_column(3) / 1000000 / 2
# Qe_k = kobra.get_column(2) / 1000000 * np.sqrt(18.36)
# Qi_k = kobra.get_column(4) / 1000000 * np.sqrt(18.36)
# fe_k = Qe_k
# fi_k += Qi_k * 1.07

import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

# 设置全局字体
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["mathtext.fontset"] = "cm"  # 数学字体和 Times 更协调
plt.rcParams["font.size"] = 25
plt.rcParams["font.weight"] = "bold"  # 全局加粗

# 创建图像和坐标轴
fig, ax1 = plt.subplots(figsize=(7, 5))

labels = [
    r"$e^{-}$ Kobra",
    r"$ion^{+}$ Kobra",
    r"$ion^{+}$ FSM",
    r"$ion^{+}$ VPM",
]

# 绘图
ax1.plot(t_k, fe_k, "r-", linewidth=1.8, label=labels[0])
ax1.plot(t_k, fi_k, "b-.", linewidth=1.8, label=labels[1])
ax1.plot(t_fsm, f_fsm, "g-.", linewidth=1.8, label=labels[2])
ax1.plot(t_vp, f_vp, "k--", linewidth=1.8, label=labels[3])

# 坐标轴标签和范围
ax1.set_ylabel(r"$\Gamma_{\parallel}$ ($n_{ped}v_{Te}$)", fontweight="bold")
ax1.set_xlabel(r"t($\tau_e$)", fontweight="bold")
ax1.set_xlim(0, 4)
ax1.set_ylim(0, 0.025)

# 设置 x 轴主刻度为 1, 2, 3, 4
ax1.set_xticks([0, 1, 2, 3, 4])
# 设置 y 轴主刻度为 1, 2, 3, 4
# ax1.set_yticks([0.00, 0.02, 0.04, 0.06, 0.08])

# 添加次刻度（在两个主刻度之间再加一个小刻度）
ax1.xaxis.set_minor_locator(AutoMinorLocator(2))
ax1.yaxis.set_minor_locator(AutoMinorLocator(2))

# 去掉网格
ax1.grid(False)

# 显示四周边框 + 加粗
for spine in ax1.spines.values():
    spine.set_linewidth(1.5)

# 主刻度线样式
ax1.tick_params(
    which="major", direction="in", length=6, width=1.2, top=True, right=True
)

# 次刻度线样式
ax1.tick_params(
    which="minor", direction="in", length=3, width=1.2, top=True, right=True
)

# 图例（加粗 + Times New Roman）
ax1.legend(
    loc="upper left",
    frameon=False,
    prop={"family": "Times New Roman", "size": 20, "weight": "bold"},
)

# 标题（可选）
plt.title("(A)", fontweight="bold")

# 布局 & 保存
plt.tight_layout()
plt.savefig("ParticleFlux.pdf", bbox_inches="tight", transparent=False)  # 保存为 PDF
plt.show()
