import Freestreamingmodel as FS
import KobraOutput as KOP
import VPcode as VP
import numpy as np
import matplotlib.pyplot as plt


fsm = FS.FSM(
    FS.Setting(True, False, True, False), "Human-defined Model", FS.Ion("H+", 18.36, 1)
)
kobra = KOP.Kobra("/Users/wangce/Vlasov/VPcode/plot/flux.out")
vpcode = VP.VPcode("/Users/wangce/Vlasov/VPcode/runtime_diags.csv", 10)

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
Qe_k = kobra.get_column(1) / 1000000 / 2
Qi_k = kobra.get_column(4) / 1000000 * np.sqrt(18.36)
fe_k =  Qe_k * 0.9
fi_k += Qi_k * 1.07

# 创建图像和坐标轴
fig, ax1 = plt.subplots(figsize=(10, 5))
labels = [
    r"$e^{-}$ flux on Kobra",
    r"$H^{+}(m=18.36)$ flux on Kobra",
    r"Analytical $H^{+}(m=18.36)$ flux",
    r"VPcode $H^{+}(m=18.36)$ flux",
]

# 左轴绘图(Kobra)
ax1.plot(t_k, fe_k, "r-", label=labels[0])
ax1.plot(t_k, fi_k, "b-", label=labels[1])
ax1.set_ylabel(r"Kobra Flux")
ax1.set_xlabel(r"Time/$\tau_e$")
ax1.set_xlim(0, 0.6)
ax1.set_ylim(0, 0.0005)
ax1.tick_params(axis="y", labelcolor="black")
ax1.grid(True)

# 右轴绘图(FS VP)
ax2 = ax1.twinx()
ax2.plot(t_fsm, f_fsm, "b-.", label=labels[2])
ax2.plot(t_vp, f_vp, "b--", label=labels[3])
ax2.set_ylim(0, 0.0005)
ax2.set_ylabel(r"Analytical Flux($n_ec_s$)")
ax2.tick_params(axis="y", labelcolor="black")

# 合并图例
lines_1, labels_1 = ax1.get_legend_handles_labels()
lines_2, labels_2 = ax2.get_legend_handles_labels()
ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc="upper left")

draw = True
# 显示与保存
plt.title(" ")
plt.tight_layout()
plt.savefig("Comparison.png", dpi=300, transparent=True)
plt.show()
