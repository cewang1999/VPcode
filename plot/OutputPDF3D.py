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
ax = fig.add_subplot(111, projection='3d')

# 设置 z 的范围
z_min, z_max = np.min(data), np.max(data)

# 初始设置图形属性
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
# ax.set_zlim(z_min, z_max)

# 设置从俯视角度看
ax.view_init(elev=90, azim=0)  # elev=90 表示从上方俯视

# 更新函数
def update(frame):
    ax.clear()  # 清空当前绘制内容
    # ax.set_zlim(z_min, z_max)  # 重新设置 z 轴范围
    ax.plot_surface(X, Y, data[frame], cmap='plasma')

    # 提取文件名中的编号
    fname = os.path.basename(vlasov_files[frame])
    step = fname.split('_')[1].split('.')[0]  # 提取时间步编号
    ax.set_title(f"Time step: {int(step)}")

# 创建动画
ani = FuncAnimation(fig, update, frames=len(vlasov_files), interval=100)

plt.show()
