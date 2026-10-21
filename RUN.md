# 工作站运行说明

## 编译与运行
```bash
make                      # g++ -O3 -fopenmp
./run_case.sh 4 16        # Case4，16 线程，后台运行，结果在 runs/Case4/
tail -f runs/Case4/log.txt
```
也可以直接运行 `./vp <case编号>`（在含有 `case/`、`input/`、`results/`、`Grid/` 的目录下）。

## 对比算例（λ = 1e-3，Nx = 1000 即 dx = 2λ，Nv = 256，0–4 τi）

| Case | dt | x-CFL(v=6) | v 子循环 | 壁面出射面值 | 用途 |
|---|---|---|---|---|---|
| 4 | 0.004 | 12 | 无 | PFC3 重构 | 文献分辨率的 AP 基准（当前格式） |
| 5 | 2.5e-4 | 0.75 | VCFL 0.9 | upwind1 = 贴壁网格值 | 仿 kobra：x、v 方向 CFL < 1，壁面一阶 |
| 6 | 2.5e-4 | 0.75 | VCFL 0.9 | PFC3 重构 | 与 Case5 只差壁面处理 |

在本机 (6 线程) 测得 Case5 约 0.06 s/步，4 τi 共 685600 步约 11 h；Case4 约 42850 步。
只看上升段和峰值可以把 Case5/6 的 `num_ts` 改为 205700（1.2 τi）。

## case 文件可选参数（缺省值）
- `lambda`（1e-3）：λ_D / L_par
- `WallFace`（PFC3）：`PFC3` 或 `upwind1`（要求 x-CFL < 1）
- `VCFL`（0）：v 方向子循环的 CFL 上限，0 为不子循环
- `OutDiag`（10）：`runtime_diags.csv` 输出间隔（步）
- `OutF`（100）：电子分布函数输出间隔（步），0 为不输出

## runtime_diags.csv 各列（每个粒子种类一组）
- `ts`：时间 (τi)
- `flux.<sp>` / `Enflux.<sp>`：壁面净粒子流 / 能流（n v_Ti，n v_Ti T），两壁平均
- `flux1.<sp>`：与净通量同一时刻的“贴壁网格值 × v”（出射方向），与 kobra 的 wallinfo 相同
- `N.<sp>`：域内总粒子数（n L），可检查 N(0) − N(t) = ∫ 2Γ dt 以及 N_i − N_e
- `den.<sp>`：贴壁网格密度；`totalKE.<sp>`：总动能；最后一列 `ElectricEnergy`

画图：`~/plotenv/bin/python plot/FluxComparison.py runs/Case5/runtime_diags.csv runs/Case5`
