import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys


class VPcode:
    def __init__(self, filepath, Col):
        self.filepath = filepath
        self.Col = Col

    def readFiles(self):
        # 如果提供了命令行参数
        if len(sys.argv) > 1:
            self.filepath = sys.argv[1]
        if len(sys.argv) > 2:
            self.Col = int(sys.argv[2])
        # === 读取CSV数据 ===
        data = pd.read_csv(self.filepath)
        # 检查列数是否足够
        if self.Col >= data.shape[1]:
            print(f"Error: CSV 只有 {data.shape[1]} 列，无法读取第 {self.Col} 列")
            sys.exit(1)
        return data

    def getTime(self):
        data = self.readFiles()
        return data.iloc[:, 0]

    def getCol(self):
        data = self.readFiles()
        return data.iloc[:, self.Col]

    def getColRequire(self, num):
        data = self.readFiles()
        return data.iloc[:, num]

    def draw(self):
        time = self.getTime()
        value = self.getCol()

        plt.figure(figsize=(10, 5))
        plt.plot(time, value, label=f"Column {self.Col}")
        plt.xlabel("Time")
        plt.ylabel("Value")
        # plt.xlim(0, 4)
        # plt.ylim(0, 0.004)
        plt.title(f"Time vs Column {self.Col}")
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.show()


# === 参数设置 ===
if __name__ == "__main__":
    csv_file = (
        "/Users/wangce/Vlasov/remote/HighResOutput/runtime_diags.csv"  # 默认CSV文件名，可用命令行替换
    )
    column_index = 2  # 要绘制的列索引（从1开始，因为0是时间）

    vpcode = VPcode(csv_file, column_index)
    vpcode.draw()