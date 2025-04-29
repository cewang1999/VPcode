import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# --- 设置路径
grid_file = os.path.join('Grid', 'Grid.txt')
vlasov_files = sorted(glob.glob(os.path.join('results', 'vlasov_*.txt')))

# --- 读取网格
grid = []
with open(grid_file, 'r') as f:
    for line in f:
        row = list(map(float, line.strip().split(',')))
        grid.append(row)

# x, y 数据应该分别为每列和每行的唯一值
x = np.array(grid[0])  # 第一行
y = np.array(grid[1])  # 第二行

# 推断行列数
unique_x = np.unique(x)
unique_y = np.unique(y)

cols = len(unique_x)
rows = len(unique_y)

print(f"推断出的网格尺寸: rows={rows}, cols={cols}")

# 创建网格
X, Y = np.meshgrid(unique_x, unique_y)  # 使用唯一的 x 和 y 创建网格

# --- 读取数据
frames = []
for fname in vlasov_files:
    z = np.loadtxt(fname, delimiter=',')
    frames.append(z.reshape((rows, cols)))

data = np.stack(frames, axis=0)  # (time_steps, rows, cols)

# --- 动画绘制
fig = plt.figure()
ax = fig.add_subplot(111)

# 设置等高线图的范围
z_min, z_max = np.min(data), np.max(data)

# 初始设置图形属性
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_title('Contour Plot')

# 只创建一次颜色条
CS = ax.contour(X, Y, data[0], cmap='viridis')
cbar = fig.colorbar(CS)

# 更新函数
def update(frame):
    ax.clear()  # 清空当前绘制内容
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_title(f"Time step: {frame}")
    
    # 绘制等高线
    CS = ax.contour(X, Y, data[frame], cmap='plasma')
    
    # 更新颜色条
    # cbar.set_mappable(CS.collections[0])  # 重新设置颜色条的 mappable

# 创建动画
ani = FuncAnimation(fig, update, frames=len(vlasov_files), interval=100)

plt.show()
