import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys

# === 参数设置 ===
csv_file = "runtime_diags.csv"  # 默认CSV文件名，可用命令行替换
column_index = 9  # 要绘制的列索引（从1开始，因为0是时间）

# 如果提供了命令行参数
if len(sys.argv) > 1:
    csv_file = sys.argv[1]
if len(sys.argv) > 2:
    column_index = int(sys.argv[2])

# === 读取CSV数据 ===
data = pd.read_csv(csv_file)

# 检查列数是否足够
if column_index >= data.shape[1]:
    print(f"Error: CSV 只有 {data.shape[1]} 列，无法读取第 {column_index} 列")
    sys.exit(1)

# === 提取时间和所选列 ===
time = data.iloc[:, 0] 
value = data.iloc[:, column_index] 

# === 绘图 ===
plt.figure(figsize=(10, 5))
plt.plot(time, value, label=f"Column {column_index}")
plt.xlabel("Time")
plt.ylabel("Value")
plt.title(f"Time vs Column {column_index}")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
