import matplotlib.pyplot as plt
import numpy as np

from matplotlib.ticker import AutoMinorLocator
from matplotlib.colors import LogNorm

# 设置全局字体
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["mathtext.fontset"] = "cm"  # 数学字体和 Times 更协调
plt.rcParams["font.size"] = 25
plt.rcParams["font.weight"] = "bold"  # 全局加粗
plt.rcParams["mathtext.default"] = "regular"  # LaTeX 数学模式使用正体

# 生成示例二维数据
x = np.linspace(-1, 1, 100)
y = np.linspace(-8, 8, 100)
X, Y = np.meshgrid(x, y)
Z = np.exp(-(X**2 + Y**2))  # 高斯函数示例

# 归一化到0~1
Z = Z / Z.max()
Z = Z * 4

# 创建图形和坐标轴
fig, ax = plt.subplots(figsize=(7, 5))
fig.subplots_adjust(left=0.2, right=0.95, top=0.9, bottom=0.15)

# 绘制二维热图
c = ax.imshow(
    Z,
    extent=[x.min(), x.max(), y.min(), y.max()],
    origin="lower",
    cmap="jet",
    aspect="auto",
    vmin=1e-3,
    vmax=4.0,
)

# 坐标轴标签和范围
ax.set_title(r"$(A)(B)t=0.3\tau_e$", fontweight="bold")
ax.set_ylabel(r"$v_{x}(v_{Te})$", fontweight="bold")
ax.set_xlabel(r"x", fontweight="bold")
ax.set_xlim(-1, 1)
ax.set_ylim(-8, 8)

# 设置 x 轴主刻度为 1, 2, 3, 4
ax.set_xticks([-1.0, -0.5, 0, 0.5, 1])
# 设置 y 轴主刻度为 1, 2, 3, 4
ax.set_yticks([-8, -6, -4, -2, 0, 2, 4, 6, 8])


# 添加次刻度（在两个主刻度之间再加一个小刻度）
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.yaxis.set_minor_locator(AutoMinorLocator(5))

# 显示四周边框 + 加粗
for spine in ax.spines.values():
    spine.set_linewidth(1.5)

# 主刻度线样式
ax.tick_params(
    which="major", direction="out", length=6, width=1.2, top=True, right=True
)

# 次刻度线样式
ax.tick_params(
    which="minor", direction="out", length=3, width=1.2, top=True, right=True
)


# 自定义 colorbar 显示范围
# vmin, vmax = 1, 4  # 只显示这个范围的颜色
color_ticks = [1, 2, 3, 4]  # 自定义刻度
cb = fig.colorbar(c, ax=ax, ticks=color_ticks)
# cb.set_label("CB")

# 布局 & 保存
# plt.tight_layout()
# plt.savefig(
#     "FastElectronDistri.pdf", bbox_inches="tight", transparent=False
# )  # 保存为 PDF
plt.show()
