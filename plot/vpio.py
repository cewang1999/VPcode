"""VPcode 输出读取 (只依赖 numpy)。

VPD 格式 (VPcode dump, 第 1 版)
--------------------------------
文件 = 文本头 + 二进制数据。文本头为 UTF-8 文本行，可直接 `head -c 4000 文件` 查看：

    VPD 1
    meta <键> <值...>            # 元数据: case、git、dt、λ、单位、case 文件参数 (param.<键>) 等
    fix  <名> <类型> <维度...>    # 只存一次的数组 (网格)
    rec  <名> <类型> <维度...>    # 每条记录都有的数组 (标量无维度)
    end

`end\\n` 之后依次为全部 fix 数组，然后是若干条记录，每条依次含全部 rec 数组。
C 序，字节序见 meta endian，类型 f8 / f4 / i8。记录数 = (文件长 − 头 − fix) / 单条长度，
末尾不完整的记录 (算例中断) 被忽略。

一个算例目录 (runs/CaseN) 下的文件:
    out/f_<ts>.vpd    快照 (1 条记录): ts, t, phi, E, rho, f_e[ix,iv], f_i[ix,iv]
    out/prof.vpd      剖面时间序列: ts, t, phi, E, rho, n/u/T/q_<e|i>, EE
    runtime_diags.csv 全局诊断 (壁面通量、能流、壁面 n/T 等)
旧算例 (results/vlasov_*.txt + Grid/Grid.txt，仅电子) 由 Run 自动兼容读取。

用法:
    from vpio import Run
    r = Run("runs/Case63")
    p = r.prof                      # p["t"] (nrec,), p["T_e"] (nrec, nx), p["x"] (nx,)
    s = r.snap(t=1.0)               # 离 t=1 τi 最近的快照: s["f_e"] (nx, nv), s["v_e"], s["E"] ...
    d = r.diag                      # d["t"], d["flux.e-"], d["Tw.H+"] ...
命令行: python vpio.py runs/Case63 [文件.vpd]   打印算例 / 文件概要
"""

import glob
import os
import re
import sys

import numpy as np

_TYPES = {"f8": "f8", "f4": "f4", "i8": "i8"}


def _num(s):
    for conv in (int, float):
        try:
            return conv(s)
        except ValueError:
            pass
    return s


class VPD(dict):
    """一个 VPD 文件：dict 形式访问数组；.meta 元数据，.param case 参数，.nrec 记录数。

    rec 数组的第一维为记录 (nrec, ...)；squeeze=True 且只有 1 条记录时去掉这一维 (快照)。
    大数组以只读内存映射读取，用到的部分才真正读盘。
    """

    def __init__(self, path, squeeze=None, mmap=True):
        super().__init__()
        self.path = path
        with open(path, "rb") as fh:
            head = b""
            while b"\nend\n" not in head:
                chunk = fh.read(1 << 16)
                if not chunk:
                    raise ValueError(f"{path}: 不是完整的 VPD 文件 (找不到 end)")
                head += chunk
        hlen = head.index(b"\nend\n") + 5
        lines = head[:hlen].decode("utf-8").splitlines()
        if not lines[0].startswith("VPD"):
            raise ValueError(f"{path}: 不是 VPD 文件")
        self.version = int(lines[0].split()[1])
        self.meta, self.param, fix, rec = {}, {}, [], []
        for ln in lines[1:]:
            parts = ln.split()
            if not parts or parts[0] == "end":
                continue
            if parts[0] == "meta":
                key, val = parts[1], ln.split(None, 2)[2] if len(parts) > 2 else ""
                if key.startswith("param."):
                    self.param[key[6:]] = _num(val)
                else:
                    self.meta[key] = _num(val)
            elif parts[0] in ("fix", "rec"):
                spec = (parts[1], _TYPES[parts[2]], tuple(int(d) for d in parts[3:]))
                (fix if parts[0] == "fix" else rec).append(spec)
        bo = "<" if self.meta.get("endian", "little") == "little" else ">"
        fdt = np.dtype([(n, bo + t, s) for n, t, s in fix])
        rdt = np.dtype([(n, bo + t, s) for n, t, s in rec])
        size = os.path.getsize(path)
        self.nrec = (size - hlen - fdt.itemsize) // rdt.itemsize if rdt.itemsize else 0
        if fdt.itemsize:
            fx = np.fromfile(path, dtype=fdt, count=1, offset=hlen)[0]
            for n, _, _ in fix:
                self[n] = np.array(fx[n])
        if self.nrec > 0:
            off = hlen + fdt.itemsize
            if mmap:
                rc = np.memmap(path, dtype=rdt, mode="r", offset=off, shape=(self.nrec,))
            else:
                rc = np.fromfile(path, dtype=rdt, count=self.nrec, offset=off)
            if squeeze is None:
                squeeze = self.meta.get("kind") == "snapshot"
            for n, _, _ in rec:
                self[n] = rc[n][0] if (squeeze and self.nrec == 1) else rc[n]
        self.fields = [n for n, _, _ in fix + rec]

    def __getattr__(self, k):
        try:
            return self[k]
        except KeyError:
            raise AttributeError(k)

    def summary(self):
        out = [f"{self.path}: VPD {self.version}, kind={self.meta.get('kind')}, 记录数 {self.nrec}"]
        for k in ("case", "git", "created", "nx", "nv", "dt", "lambda"):
            if k in self.meta:
                out.append(f"  {k} = {self.meta[k]}")
        for n in self.fields:
            a = self[n]
            out.append(f"  {n:8s} {str(a.dtype):8s} {a.shape}")
        if "t" in self and self.nrec:
            t = np.atleast_1d(self["t"])
            out.append(f"  t (τi): {t[0]:.4g} … {t[-1]:.4g}")
        return "\n".join(out)


def read_vpd(path, **kw):
    return VPD(path, **kw)


def moments(f, v, mass, axis=-1):
    """由分布函数 f(..., v) 计算 n、u、T、q (同 out/prof.vpd 的定义，v 为均匀速度格中心)。"""
    dv = v[1] - v[0]
    n = f.sum(axis) * dv
    u = (f * v).sum(axis) * dv / np.where(n > 0, n, 1)
    c = v - np.expand_dims(u, axis)
    T = mass * (f * c**2).sum(axis) * dv / np.where(n > 0, n, 1)
    q = 0.5 * mass * (f * c**3).sum(axis) * dv
    return n, u, T, q


def read_case(path):
    """读取 case/CaseN.dat 为 dict。"""
    d = {}
    for ln in open(path):
        p = ln.split()
        if p and not p[0].startswith("#"):
            d[p[0]] = _num(p[1]) if len(p) > 1 else ""
    return d


def read_diag(path):
    """读取 runtime_diags.csv 为 dict (列名 → 数组)；第 0 列 (时间, τi) 另存为 "t"。"""
    with open(path) as fh:
        names = fh.readline().strip().split(",")
    a = np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2)
    d = {n: a[:, k] for k, n in enumerate(names)}
    d["t"] = a[:, 0]
    return d


class Run:
    """一个算例目录。属性均为惰性读取。"""

    MASS_RATIO = 1836.0

    def __init__(self, path):
        self.path = path.rstrip("/")
        self.out = os.path.join(self.path, "out")
        self._prof = self._diag = None
        self.snap_files = sorted(glob.glob(os.path.join(self.out, "f_*.vpd")))
        self.legacy = not self.snap_files and bool(glob.glob(os.path.join(self.path, "results", "vlasov_*.txt")))
        cases = glob.glob(os.path.join(self.path, "case", "Case*.dat"))
        self.case = read_case(cases[0]) if cases else {}
        if self.legacy:
            self.snap_files = sorted(glob.glob(os.path.join(self.path, "results", "vlasov_*.txt")))

    # ---- 剖面与诊断 ----
    @property
    def prof(self):
        if self._prof is None:
            p = os.path.join(self.out, "prof.vpd")
            if not os.path.exists(p):
                raise FileNotFoundError(f"{p} 不存在 (旧算例或 OutProf 0)")
            self._prof = VPD(p, squeeze=False)
        return self._prof

    @property
    def diag(self):
        if self._diag is None:
            self._diag = read_diag(os.path.join(self.path, "runtime_diags.csv"))
        return self._diag

    # ---- 快照 ----
    def _ts(self, fn):
        return int(re.findall(r"(\d+)\.(?:vpd|txt)$", fn)[0])

    @property
    def snap_ts(self):
        return np.array([self._ts(f) for f in self.snap_files], dtype=int)

    @property
    def snap_t(self):
        """各快照的时间 (τi)。"""
        if self.legacy:
            # 旧格式 vlasov_<ts>.txt 为第 ts 步推进之后的 f，即 f^{ts+1}
            return (self.snap_ts + 1) * self.case["dt"] / np.sqrt(self.MASS_RATIO)
        return np.array([VPD(f)["t"] for f in self.snap_files])

    def snap(self, k=None, t=None):
        """第 k 个快照，或时间离 t (τi) 最近的快照。"""
        if t is not None:
            k = int(np.argmin(np.abs(self.snap_t - t)))
        fn = self.snap_files[k]
        if not self.legacy:
            return VPD(fn)
        return self._legacy_snap(fn)

    def _legacy_snap(self, fn):
        g = os.path.join(self.path, "Grid", "Grid.txt")
        with open(g) as fh:
            x = np.array(fh.readline().split(","), float)
            v = np.array(fh.readline().split(","), float)
        s = VPD.__new__(VPD)
        dict.__init__(s)
        s.path, s.version, s.nrec, s.param = fn, 0, 1, self.case
        s.meta = {"kind": "legacy", "nx": len(x), "nv": len(v), "dt": self.case.get("dt")}
        s["x"], s["v_e"] = x, v
        s["ts"] = self._ts(fn) + 1
        s["t"] = s["ts"] * self.case["dt"] / np.sqrt(self.MASS_RATIO)
        s["f_e"] = np.loadtxt(fn, delimiter=",").T  # 文件为 nv 行 × nx 列 (6 位有效数字)
        s.fields = ["x", "v_e", "ts", "t", "f_e"]
        return s

    def summary(self):
        out = [f"{self.path}: " + ("旧格式 (仅电子 f, 6 位有效数字)" if self.legacy else "VPD")]
        if self.case:
            keys = ("Nx", "Nv", "dt", "SolverType", "XGrid", "WallE", "num_ts")
            out.append("  case: " + ", ".join(f"{k} {self.case[k]}" for k in keys if k in self.case))
        if self.snap_files:
            t = self.snap_t
            out.append(f"  快照 {len(t)} 个, t = {t[0]:.4g} … {t[-1]:.4g} τi")
        if os.path.exists(os.path.join(self.out, "prof.vpd")):
            p = self.prof
            out.append(f"  剖面 {p.nrec} 条, t = {p['t'][0]:.4g} … {p['t'][-1]:.4g} τi")
        if os.path.exists(os.path.join(self.path, "runtime_diags.csv")):
            d = self.diag
            out.append(f"  诊断 {len(d['t'])} 行, t 至 {d['t'][-1]:.4g} τi")
        return "\n".join(out)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(0)
    for a in sys.argv[1:]:
        print(VPD(a).summary() if a.endswith(".vpd") else Run(a).summary())
