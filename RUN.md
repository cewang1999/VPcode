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

## λ 扫描（PFC3，dx = 2λ，Δt = 4/ω_pe = 4λ，同文献；0–4 τi）

| Case | λ | Nx | dt | 步数 | 用途 |
|---|---|---|---|---|---|
| 4 | 1e-3 | 1000 | 0.004 | 42850 | 文献图 1 上排 |
| 7 | 5e-4 | 2000 | 0.002 | 85700 | 中间 λ |
| 8 | 2.5e-4 | 4000 | 0.001 | 171400 | 文献图 1 下排 |
| 9 | 2.5e-4 | 1000 (dx = 8λ) | 0.004 | 42850 | 粗网格对照 |
| 10 | 2.5e-4 | 1000 (dx = 8λ) | 0.001 (ω_pe Δt = 4) | 171400 | 推测的文献设置：1000 个 x 点 + Δt = 4/ω_pe |
| 11 | 2.5e-4 | 1000, Nv = 1024 | 0.001 | 171400 | 同 Case10，速度点数同文献 (1000) |
| 12 | 1e-3 | 1000, Nv = 1024 | 0.004 | 42850 | 同 Case4，速度点数 ×4：检查电子温度 / 图 7 |
| 13 | 1e-3 | 1000, Nv = 1024 | 0.001 (ω_pe Δt = 1) | 171400 | 同 Case12，时间步 /4：检查电子数值加热 |
| 14 | 1e-3 | 4000 (dx = 0.5λ) | 0.001 | 171400 | 加密网格：检验后期电子加热是否来自未分辨的局部德拜长度（dt=0.004 时鞘层分辨后失稳） |
| 15 | 1e-3 | 8000 (dx = 0.25λ) | 0.0005 | 171400 (到 2 τi) | 同上，再加密一倍看收敛 |
| 16 | 1e-3 | 1000 | 0.004 | 42850 | 同 Case4，`WallE onesided` |
| 17 | 1e-3 | 1000, Nv = 1024 | 0.001 | 171400 | 同 Case13，`WallE onesided` |
| 18 | 2.5e-4 | 1000, Nv = 1024 | 0.001 | 171400 | 同 Case11（加热最严重），`WallE onesided` |
| 19 | 1e-3 | 1000 | 0.004 | 42850 | 同 Case4，`WallE kobra`（壁面离散与 quasi-kobra 相同） |
| 20 | 1e-3 | 1000 | 0.004 | 42850 | 同 Case16，加 `EFilter 1`（内部电场二项式滤波）：检验后期 (t>2.8) 的栅格尺度不稳定性 |
| 21 | 1e-3 | 2000 (dx = λ) | 0.002 | 85700 | 同 Case16，网格加密一倍：不稳定性应推迟或减弱 |
| 22 | 1e-3 | 1000 | 0.008 | 21425 | 同 Case16，dt 加倍 (ω_pe Δt = 8)：AP 阻尼更强，不稳定性应减弱；若壁面速度越界会失稳 |
| 23 | 1e-3 | 1000, Nv = 1024 | 0.001 | 171400 | 同 Case17，加 `SpikeDump 5`：抓 t≈1.8–3.5 τi 电子通量 / 电场能的突发事件，结果在 `spike/` |
| 24 | 1e-3 | 1000, Nv = 256 | 0.001 | 171400 | 同 Case23，Nv = 256（快约 4 倍）：看突发事件是否与速度分辨率有关 |
| 25 | 1e-3 | 1000, Nv = 1000 | 0.004 | 42850 | 文献分辨率（x、v 各 1000 点，Δt = 4/ω_pe），`WallE onesided` + `SpikeDump 5` |
| 26 | 1e-3 | 1000, Nv = 1000 | 0.004 | 42850 | 同 Case25，`WallE central`（原壁面离散）：修复舍入错误后在文献分辨率下最终确认 |
| 27 | 1e-3 | 1000, Nv = 256 | 0.001 | 42850 | 与 kobra 1000×256 同分辨率（ω_pe Δt = 1，与 kobra 自适应步长相近），`WallE central`，0–1 τi |
| 28 | 1e-3 | 512, Nv = 256 | 0.001 | 25710 | 与 kobra 512×256 同分辨率，`WallE central`，0–0.6 τi |
| 29 | 1e-3 | 256, Nv = 128 | 0.0005 | 64275 | 与 kobra 256×128 同分辨率，`WallE central`，0–0.75 τi |
| 30 | 1e-3 | 1000, Nv = 256 | 自适应 | 0–1 τi | 与 kobra 1000×256 完全相同的离散：`Scheme kobra` + `WallE kobra` + `Umax 6` |
| 31 | 1e-3 | 512, Nv = 256 | 自适应 | 0–0.6 τi | 同 Case30，对应 kobra 512×256 |
| 32 | 1e-3 | 256, Nv = 128 | 自适应 | 0–0.75 τi | 同 Case30，对应 kobra 256×128（本机约 1 分钟） |
| 33 | 1e-3 | 1000, Nv = 256 | 自适应 | 0–1 τi | 同 Case30，`WallE central`：只换壁面离散 |
| 34 | 1e-3 | 1000, Nv = 256 | 自适应 | 0–1 τi | 同 Case30，`Flux PFC3`：只换输运格式（本机约 3 h，建议工作站） |
| 35 | 1e-3 | 1000, Nv = 256 | 自适应 | 0–1 τi | 同 Case30，`Flux PFC3` + `WallE central`：与 Case27 只差时间推进（建议工作站） |
| 36–38 | 1e-3 | 512, Nv = 256 | 自适应 | 0–0.6 τi | 同 Case33–35，对应 Case31（本机约 40 min） |
| 39 | 1e-3 | 512, Nv = 256 | 自适应 | 0–0.6 τi | 同 Case37（PFC3 + WallE kobra），`WallFace upwind1`：贴壁格偏高 +34% → +6.5%，与 kobra 相同，原因在壁面出射面值 |
| 40 | 1e-3 | 512, Nv = 256 | 自适应 | 0–0.6 τi | 同 Case36（MUSCL），`WallE onesided` |
| 41 | 1e-3 | 1000, Nv = 256 | 自适应 | 0–1 τi | 同 Case33（MUSCL），`WallE onesided` |
| 42 | 1e-3 | 512, Nv = 256 | 自适应 | 0–0.6 τi | 同 Case38（PFC3），`WallE onesided`：onesided 下后期升温消失，与输运格式无关 |

## case 文件可选参数（缺省值）
- `lambda`（1e-3）：λ_D / L_par
- `WallFace`（PFC3）：`PFC3` 或 `upwind1`（要求 x-CFL < 1）
- `VCFL`（0）：v 方向子循环的 CFL 上限，0 为不子循环
- `OutDiag`（10）：`runtime_diags.csv` 输出间隔（步）
- `OutF`（100）：电子分布函数输出间隔（步），0 为不输出
- `EFilter`（0）：内部电场 (1,2,1)/4 二项式滤波次数，贴壁 2 格不滤
- `SpikeDump`（0）：电场能突跳诊断。某步电场能超过上一步的 `SpikeDump` 倍时，在 `spike/ev<k>_` 下保存触发步的 fⁿ 与前一步 f^{n-1}（`fe_now/fe_prev/fi_now/fi_prev.txt`，格式同 `results/vlasov_*.txt`），以及触发前后各 20 步的场剖面 `prof.csv`（φ, E, ρ, J, R, Σf/m，含 i = −1 与 Nx 两个幽灵点）；最多 5 次
- `WallE`（central）：贴壁格电场离散。`central`：原方案，E₀ = −(φ₁ − φ_g)/2dx，φ_g = −φ₀；`onesided`：E₀ = −(φ₁ − φ₀)/dx，且 AP 方程壁面面系数只保留 λ²（满足离散高斯定律，后期电子数值加热大幅减小）；`kobra`：与 quasi-kobra 的 mug_kobra.c 相同（幽灵 φ_g = −2.5φ₀ + 0.5φ₁，壁面电流取壁面面通量，壁面 ∂R 单侧估计）
- `Scheme`（PFC3）：输运与时间推进。`kobra`：与 quasi-kobra 相同——x、v 方向均为 MUSCL-Hancock（kobra 的 van Albada 斜率 `av`，即 `semi`），
  x 吸收壁按 kobra 的 `usrbd`（出射面值 = 贴壁格值，入射面值 = 下限），v 方向 FREE 边界；
  每步 dt = `CFL` / max(|v|/dx, |qE/m|/dv)（只统计 f 高于下限的格，同 `usrtimestep`），由 fⁿ 解 φ 后做 Strang 分裂 x(dt/2) → v(dt) → x(dt/2)；
  f 不截断而是抬到下限 `Dfloor`·m·λ²（同 `usrfxp`）。算到 t_end = `dt` × `num_ts`（case 中的 dt 只用来定 t_end）。
  此时 `clipNeg/clipSmall` 记录的是被抬到下限所加的质量，`flux1/nw/Tw` 取每步结束时贴壁格出射部分（同 kobra 的诊断）
- `Flux`（MUSCL）：`Scheme kobra` 下的输运格式。`MUSCL` 同 kobra；`PFC3` 为 VPcode 原格式（壁面出射面值见 `WallFace`），
  时间推进仍同 kobra，用来把输运格式、壁面离散、时间推进三者的影响分开
- `Dfloor`（1e-10）、`CFL`（0.9）：见 `Scheme kobra`
- `Umax`（0）：速度网格 ±Umax 个热速度（电子 v_Te，离子 v_Ti），对应 kobra 的 `U … -Umax Umax`；0 为原网格（电子 ±6，离子 ±0.15 v_Te = ±6.43 v_Ti）

## runtime_diags.csv 各列（每个粒子种类一组）
- `ts`：时间 (τi)
- `flux.<sp>` / `Enflux.<sp>`：壁面净粒子流 / 能流（n v_Ti，n v_Ti T），两壁平均
- `flux1.<sp>`：与净通量同一时刻的“贴壁网格值 × v”（出射方向），与 kobra 的 wallinfo 相同
- `N.<sp>`：域内总粒子数（n L），可检查 N(0) − N(t) = ∫ 2Γ dt 以及 N_i − N_e
- `nw.<sp>` / `Tw.<sp>`：壁面密度 / 平行温度（文献图 7），由壁面分布 f_wall = 面通量 / v（出射方向）、入射方向为 0 求得
- `den.<sp>`：贴壁网格密度（网格中心，含被反射的电子，不是壁面值）；`totalKE.<sp>`：总动能；`ElectricEnergy` = ΣE²Δx（未乘 λ²/2）
- `clipNeg.<sp>` / `clipSmall.<sp>`（加在最后，原有列号不变）：累计被 `f < 1e-10 → 0` 截断清零的负值 / 小正值质量（n L）。PFC3 保正，`clipNeg` 应约为 0（~1e-25）；明显非零说明通量出错（如修复前的出发格舍入错误，离子 `clipNeg` 一步跳到 1e-5）

画图：
- `~/plotenv/bin/python plot/PaperFigs.py runs/Case4 [runs/Case12 ...]`：按文献坐标轴画图 1、图 7 和 Γe/Γi − 1，图存在运行目录（多个目录时存在 runs/）
- `~/plotenv/bin/python plot/FluxComparison.py runs/Case5/runtime_diags.csv runs/Case5`

## 与 kobra 对比
quasi-kobra-sup 的 `naga/usr.c` 末尾加了 `diag_mon()`（`DIAG BEGIN` 到 `DIAG END`，另在 `usrmon()` 里有一行标 `DIAG` 的调用），
每步输出 `$kobra_data/diag_budget.dat`（N、平行动能、W、高斯定律 Q 与 E_R−E_L、壁面贴壁格出射 n/j/T、floor 量、Poisson 残差），
电场能突跳时写 `diag_spike_p<rank>.dat`。kobra 跑完后：
```bash
tools/fetch_kobra.sh <标签> [kobra数据目录] [kobra算例.dat]   # 拷到 runs/kobra_<标签>/
~/plotenv/bin/python plot/KobraDiag.py runs/kobra_<标签> runs/Case25   # 对比图 KobraDiag.png + 终端汇总
```
同分辨率对照：kobra 1000×256 ↔ Case27，512×256 ↔ Case28，256×128 ↔ Case29（VPcode 原格式，central）；
同离散对照（`Scheme kobra`）：↔ Case30、Case31、Case32。
VPcode 的净通量与 kobra 一样由粒子数求（−dN/dt / 2）：flux 列每 `OutDiag` 步采样一次，ω_pe Δt ~ 1 时电子通量的等离子体振荡会混叠。
图中 Γe/Γi − 1 的净通量与贴壁格分两幅画（细线原始值，粗线 0.01 τi 滑动平均），时间轴只到 kobra 数据结束处。
kobra 的净通量由 Γ = −dN/dt / 2 得到（内部守恒，两壁平均）；kobra 为单精度，离子到达壁面前的离子"净通量"是舍入噪声。
