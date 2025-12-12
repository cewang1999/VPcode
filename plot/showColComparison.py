import KobraOutput as KOP
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

nu0 = KOP.Kobra("/Users/wangce/Vlasov/remote/ColComparsion/nu0.out")
nu2 = KOP.Kobra("/Users/wangce/Vlasov/remote/ColComparsion/nu2.out")
nu1 = KOP.Kobra("/Users/wangce/Vlasov/remote/ColComparsion/nu1.out")
Kbr = KOP.Kobra("/Users/wangce/Vlasov/remote/flux3.out")

# Kobra

t = Kbr.get_column(0) / 42.848
je = Kbr.get_column(1) / 2 * 42.848
ji = Kbr.get_column(3) / 2 * 42.848

t0 = nu0.get_column(0)
t1 = nu1.get_column(0)
t2 = nu2.get_column(0)
je0 = nu0.get_column(4)
je2 = nu2.get_column(4)
je1 = nu1.get_column(4)

# 设置全局字体
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["mathtext.fontset"] = "cm"  # 数学字体和 Times 更协调
plt.rcParams["font.size"] = 25
plt.rcParams["font.weight"] = "bold"  # 全局加粗
plt.rcParams["mathtext.default"] = "regular"  # LaTeX 数学模式使用正体

# 创建图像和坐标轴
fig, ax1 = plt.subplots(figsize=(7, 5), constrained_layout=True)

labels = [
    r"$\nu=0.0$",
    r"$\nu=0.2$",
    r"$\nu=0.1$",
]

# 绘图
ax1.plot(t, je, "r-", linewidth=1.8, label=labels[0])
ax1.plot(t, ji, "b-", linewidth=1.8, label=labels[1])
# ax1.plot(t1, je1, "k-", linewidth=1.8, label=labels[2])

# 坐标轴标签和范围
ax1.set_ylabel(r"$\Gamma$", fontweight="bold")
ax1.set_xlabel(r"t($\tau_e$)", fontweight="bold")
ax1.set_xlim(0, 4)

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
# plt.title("(B)", fontweight="bold")

# # 布局 & 保存
# plt.tight_layout()
# plt.savefig(
#     "EnergyComparation.pdf", bbox_inches="tight", transparent=False
# )  # 保存为 PDF
plt.show()
