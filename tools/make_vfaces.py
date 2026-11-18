"""生成非均匀速度网格的格面表 (case 选项 VFacesE / VFacesI 读取；格式也可直接用于 kobra 的 U 网格)。

格面以本粒子热速度为单位 (电子 v_Te，离子 v_Ti，即 kobra 的 U)，关于 0 对称，范围 [-Umax, Umax]，共 N 个格。

两种模式:
  envelope  由参考算例的分布函数统计各 |U| 处需要的格宽: 对所有快照 (t ≥ tmin) 与 n > 1e-3 n_max 的格，
            在 |U - u| < 3σ 内要求 ΔU ≤ σ/k；再按需求分配 N 个格 (格宽 ∝ 需求，即 k 按 N 自动缩放)，
            格宽上限 cap (早期热分布与尾部)。输出中给出等效的 k。
  geom      解析网格: |U| < vc 内等宽 dv0，其外格宽按 ratio 等比增大 (不超过 cap)，到 Umax 为止；N 由参数决定。

用法:
  python tools/make_vfaces.py envelope runs/Case83 --sp e --n 128 --out case/vfaces/e128.txt
  python tools/make_vfaces.py envelope runs/Case83 --sp e,i --n 256 --out case/vfaces/ei256.txt   (共用网格, 同 kobra)
  python tools/make_vfaces.py geom --dv0 0.006 --vc 0.1 --ratio 1.05 --umax 6 --out case/vfaces/eg.txt
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "plot"))


def envelope(run_dir, species, k, tmin, umax, ngrid=60001):
    """各 |U| 处需要的格宽 (k = 1 时为 σ)。返回 (U, req)，U ∈ [0, umax]。"""
    from vpio import Run, moments

    r = Run(run_dir)
    U = np.linspace(0.0, umax, ngrid)
    req = np.full(U.shape, np.inf)
    for kk in range(len(r.snap_t)):
        if r.snap_t[kk] < tmin:
            continue
        s = r.snap(kk)
        for sp in species:
            m = float(s.meta.get(f"mass_{sp}", 1.0 if sp == "e" else 1836.0))
            f = np.asarray(s[f"f_{sp}"], float)
            v = np.asarray(s[f"v_{sp}"])
            dv = np.diff(np.asarray(s[f"vf_{sp}"])) if f"vf_{sp}" in s else None
            vv, ff = v * np.sqrt(m), f / np.sqrt(m)
            dd = None if dv is None else dv * np.sqrt(m)
            n, u, T, _ = moments(ff, vv, 1.0, dv=dd)
            sig = np.sqrt(np.maximum(T, 0.0))
            ok = n > 1e-3 * n.max()
            for uj, sj in zip(np.abs(u[ok]), sig[ok]):
                lo, hi = np.searchsorted(U, [uj - 3 * sj, uj + 3 * sj])
                req[lo:hi] = np.minimum(req[lo:hi], sj)
    return U, req / k


def faces_from_req(U, req, n, cap):
    """按需求格宽 req(U) (U ≥ 0) 分配 n 个格 (n 为偶数，对称)：格宽 ∝ req，上限 cap。返回 (faces, 缩放系数)。"""
    if n % 2:
        sys.exit("N 必须为偶数 (关于 0 对称)")
    half = n // 2

    def cells(scale):  # 格宽 = min(scale * req, cap) 时半区的格数
        w = np.minimum(scale * req, cap)
        return np.trapezoid(1.0 / w, U)

    lo, hi = 1e-6, 1e6
    for _ in range(200):  # 二分: 使半区格数 = half
        mid = np.sqrt(lo * hi)
        lo, hi = (mid, hi) if cells(mid) > half else (lo, mid)
    scale = np.sqrt(lo * hi)
    w = np.minimum(scale * req, cap)
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (1 / w[1:] + 1 / w[:-1]) * np.diff(U))])
    cum *= half / cum[-1]
    pos = np.interp(np.arange(half + 1), cum, U)
    pos[0], pos[-1] = 0.0, U[-1]
    return np.concatenate([-pos[::-1], pos[1:]]), scale


def geom(dv0, vc, ratio, cap, umax):
    pos = [0.0]
    while pos[-1] < umax - 1e-12:
        h = dv0 if pos[-1] < vc - 1e-12 else min(cap, (pos[-1] - pos[-2]) * ratio if len(pos) > 1 else dv0)
        pos.append(pos[-1] + h)
    pos = np.array(pos)
    pos *= umax / pos[-1]  # 末格拉伸到正好 Umax
    return np.concatenate([-pos[::-1], pos[1:]])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mode", choices=["envelope", "geom"])
    ap.add_argument("run", nargs="?", help="envelope: 参考算例目录 (含 out/f_*.vpd)")
    ap.add_argument("--sp", default="e", help="envelope: 按哪些粒子的需求 (e、i 或 e,i)")
    ap.add_argument("--n", type=int, default=128, help="envelope: 总格数 (偶数)")
    ap.add_argument("--tmin", type=float, default=0.15, help="envelope: 只统计 t ≥ tmin (τi) 的快照")
    ap.add_argument("--cap", type=float, default=0.25, help="格宽上限 (热速度单位)")
    ap.add_argument("--umax", type=float, default=6.0)
    ap.add_argument("--dv0", type=float, default=0.006, help="geom: 中心格宽")
    ap.add_argument("--vc", type=float, default=0.1, help="geom: 中心等宽区半宽")
    ap.add_argument("--ratio", type=float, default=1.05, help="geom: 等宽区外相邻格宽之比")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    if a.mode == "envelope":
        if not a.run:
            sys.exit("envelope 模式需要参考算例目录")
        U, req = envelope(a.run, a.sp.split(","), 1.0, a.tmin, a.umax)
        req = np.where(np.isfinite(req), req, np.inf)
        faces, scale = faces_from_req(U, req, a.n, a.cap)
        note = (f"envelope from {a.run}, species {a.sp}, t >= {a.tmin} tau_i: dU = min({scale:.4g} * sigma_min(U), {a.cap}); "
                f"equivalent criterion dU <= sigma/{1 / scale:.3g}")
    else:
        faces = geom(a.dv0, a.vc, a.ratio, a.cap, a.umax)
        note = f"geom: dU0 {a.dv0} for |U| < {a.vc}, ratio {a.ratio}, cap {a.cap}, Umax {a.umax}"

    w = np.diff(faces)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w") as fh:
        fh.write(f"# 速度格面 (热速度单位，同 kobra 的 U)，N = {len(w)}；{note}\n")
        fh.write(f"# 格宽: 最小 {w.min():.5g} (|U|≈0)，最大 {w.max():.5g}\n")
        for x in faces:
            fh.write(f"{x:.12g}\n")
    print(f"{a.out}: N = {len(w)}, dU {w.min():.4g} – {w.max():.4g}；{note}")
    for u0 in (0.0, 0.1, 0.3, 0.5, 1.0, 2.0, 3.0, 5.0):
        j = np.searchsorted(faces, u0) - 1
        j = min(max(j, 0), len(w) - 1)
        print(f"   |U| = {u0:4.1f}: dU = {w[j]:.4g}")


if __name__ == "__main__":
    main()
