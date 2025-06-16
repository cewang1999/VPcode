import numpy as np


class Kobra:
    def __init__(self, filepath, skip_header=0):
        self.filepath = filepath
        self.skip_header = skip_header  # 默认不跳过任何行

    def read_out_file(self):
        """
        读取 .out 文件为 numpy 数组
        :return: numpy.ndarray，数据数组
        """
        try:
            data = np.loadtxt(self.filepath, skiprows=self.skip_header)
            return data
        except Exception as e:
            print(f"读取文件失败: {e}")
            return None

    def get_column(self, col_index):
        """
        提取指定列的 numpy 数组
        :param col_index: 列索引，从0开始
        :return: numpy.ndarray 或 None
        """
        data = self.read_out_file()
        if data is not None:
            if 0 <= col_index < data.shape[1]:
                return data[:, col_index]
            else:
                raise IndexError(
                    f"列索引 {col_index} 超出范围，最大为 {data.shape[1] - 1}"
                )
        return None


if __name__ == "__main__":
    k = Kobra("/Users/wangce/Vlasov/VPcode/plot/flux.out")  # 第一行就是数据，无需跳过
    time_data = k.get_column(0)  # 获取时间列
    exp_data = k.get_column(1)  # 获取实验数据列
    print("时间列前5个值:", time_data[:5])
    print("实验列前5个值:", exp_data[:5])
