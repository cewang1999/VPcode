import Freestreamingmodel as FS
import KobraOutput as KOP
import numpy as np
import matplotlib.pyplot as plt


fsm = FS.FSM(FS.Setting(True, False), "Human-defined Model", FS.Ion("H+", 100, 1))
kobra = KOP.Kobra("/Users/wangce/Vlasov/VPcode/plot/flux.out")

# FSM
time = np.linspace(0.01, 40, 1000)
f_fsm = fsm.Output(time)
t_fsm = fsm.TimeTransfer(time)

# Kobra
t_k = kobra.get_column(0)
fe_k = kobra.get_column(1)/10000
fi_k = kobra.get_column(3)/10000

# 创建图像和坐标轴
fig, ax1 = plt.subplots(figsize=(10, 5))
labels = [
    r"$e^{-}$ flux on Kobra",
    r"$H^{+}(m=100)$ flux on Kobra",
    r"Analytical $H^{+}(m=100)$ flux",
]

# 左轴绘图
ax1.plot(t_k, fe_k, "b-", label=labels[0])
ax1.plot(t_k, fi_k, "g--", label=labels[1])
ax1.set_ylabel(r"Kobra Flux")
ax1.set_xlabel(r"Time/$\tau_e$")
ax1.set_xlim(0, 4)
ax1.tick_params(axis="y", labelcolor="black")
ax1.grid(True)

# 右轴绘图
ax2 = ax1.twinx()
ax2.plot(t_fsm, f_fsm, "r-.", label=labels[2])
ax2.set_ylabel(r"Analytical Flux($n_ec_s$)")
ax2.tick_params(axis="y", labelcolor="black")

# 合并图例
lines_1, labels_1 = ax1.get_legend_handles_labels()
lines_2, labels_2 = ax2.get_legend_handles_labels()
ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc="upper left")

# 显示与保存
plt.title(" ")
plt.tight_layout()
plt.savefig("Comparison.png", dpi=300, transparent=True)
plt.show()
