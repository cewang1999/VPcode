# 工作站运行说明

## 编译与运行
```bash
make                      # g++ -O3 -fopenmp
./run_case.sh 4 16        # Case4，16 线程，后台运行，结果在 runs/Case4/
tail -f runs/Case4/log.txt
```
也可以直接运行 `./vp <case编号>`（在含有 `case/`、`input/` 的目录下，输出写到 `out/` 与 `runtime_diags.csv`）。

## 输出文件（VPD 格式，读取见 `plot/vpio.py`）
- `out/f_<ts>.vpd`（每 `OutF` 步）：电子、离子分布函数 f[ix][iv] 与 φ、E、ρ 剖面，以及网格（格面、格中心、格宽、两种粒子的 v）
- `out/prof.vpd`（每 `OutProf` 步追加一条）：φ、E、ρ 与电子、离子的 n、u、T、q 剖面及电场能，用于细时间分辨的剖面演化（不必输出 f）
- 文件为文本头 + 二进制：`head -c 4000 out/prof.vpd` 可直接看到 case 参数、git 版本、单位与数组列表；算例中断时已写的记录完整可读
- 写出时刻为每步解出场之后、输运之前：f 为 fⁿ（t = ts·Δt），φ/E/ρ 由同一 fⁿ 求得。注意 `runtime_diags.csv` 的 `ElectricEnergy` 在解场之前写出，是上一步的电场
- Python：`from vpio import Run; r = Run("runs/Case80"); r.prof["T_e"]`（nrec × nx）、`r.snap(t=1.0)["f_i"]`、`r.diag["flux.e-"]`；
  `python plot/vpio.py runs/Case80` 打印概要。旧算例（旧版代码输出的 `results/vlasov_*.txt` + `Grid/Grid.txt`，如 Case50、Case63）由 `Run` 自动兼容读取（仅电子 f，旧文本只有 6 位有效数字，vlasov_<ts> 为第 ts 步之后的 f^{ts+1}）

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
| 43 | 1e-3 | 1000, Nv = 1000 | 自适应 | 0–4 τi | 与 kobra 1000×1000 完全相同的离散（同 Case30，Nv = 1000），工作站运行 |
| 44 | 1e-3 | 1000, Nv = 1000 | 0.004 | 42850 | 文献方法：PFC3 + `WallE central` + 原时间推进（先解 φ^{n+1}，再 v、x），Δt = 4/ω_pe；同 Case26，只加 `Umax 6` 使速度网格与 kobra 相同，工作站运行 |
| 45 | 1e-3 | 1000, Nv = 1000 | 0.004 | 42850 | 同 Case44，加 `WallMoment out`：只把壁面面系数 / 面电流 / ∂R 换成出射方向面值（同 kobra），φ 幽灵格与 E₀ 仍为 central，工作站运行 |
| 46 | 1e-3 | 4000, Nv = 256 | 0.0005 | 342800 | 加密网格分辨鞘层（Δx = 0.5λ）：文献方法（PFC3 + `WallE central` + 原时间推进）+ `Umax 6`，0–4 τi，工作站约 3 h |
| 47 | 1e-3 | 8000, Nv = 256 | 0.0005 | 342800 | 同 Case46，Δx = 0.25λ，工作站约 6 h |
| 48 | 1e-3 | 1000, Nv = 256 | 0.0005 | 342800 | 同 Case46，Nx = 1000：同一 Δt 下的对照，三者只差 Nx |
| 49 | 1e-3 | 1000, Nv = 1000 | 0.0005 | 342800 | 同 Case48，Nv = 1000：与 Case48 只差 Nv（与 Case44 只差 Δt），工作站约 5 h |
| 50 | 1e-3 | 4000, Nv = 1000 | 0.0005 | 342800 | 同 Case46，Nv = 1000：与 Case46 只差 Nv，工作站约 6.6 h。**标杆算例**：0–4 τi 无数值升温，T_i·t² = 0.0049–0.0050，壁面通量 / 能流与 Case25 相差 < 0.3% |
| 51 | 1e-3 | 4000, Nv = 1000 | 0.002 | 85700 | Δt 扫描（以 Case50 为基准，只改 Δt）：ω_pe Δt = 2，x-CFL = 24，工作站约 1.6 h |
| 52 | 1e-3 | 4000, Nv = 1000 | 0.001 | 171400 | Δt 扫描：ω_pe Δt = 1，x-CFL = 12，约 3.3 h |
| 53 | 1e-3 | 4000, Nv = 1000 | 0.00025 | 685600 | Δt 扫描：ω_pe Δt = 0.25（显式分辨等离子体振荡），x-CFL = 3，约 13 h |
| 54 | 1e-3 | 4000, Nv = 1000 | 0.0005 | 342800 | 同 Case50，`SolverType Poisson`（标准泊松）：与 Case50 只差 AP 项，ω_pe Δt = 0.5，约 6.6 h |
| 55 | 1e-3 | 4000, Nv = 1000 | 0.001 | 171400 | 同 Case52，`SolverType Poisson`：ω_pe Δt = 1 |
| 56 | 1e-3 | 4000, Nv = 1000 | 0.002 | 85700 | 同 Case51，`SolverType Poisson`：ω_pe Δt = 2，显式稳定极限，预计发散（AP 的意义所在） |
| 57 | 1e-3 | 4000, Nv = 1000 | 0.00025 | 685600 | 同 Case53，`SolverType Poisson`：ω_pe Δt = 0.25 |
| 58 | 1e-3 | stretch → 724, Nv = 1000 | 0.0005 | 342800 | 同 Case50，非均匀 x 网格：壁面 50λ 内 0.5λ，按 1.05 渐变到体区 4λ，约 1.2 h |
| 59 | 1e-3 | stretch → 512, Nv = 1000 | 0.0005 | 342800 | 同 Case58，体区 8λ，约 0.9 h |
| 60 | 1e-3 | stretch → 422, Nv = 1000 | 0.0005 | 342800 | 同 Case58，体区 16λ（源宽 σ = 100λ 内约 6 格），约 0.7 h |
| 61 | 1e-3 | stretch → 512, Nv = 1000 | 0.001 | 171400 | 同 Case59，Δt = 0.001（Case52 表明该 Δt 在均匀网格上只使壁面 Te 低约 5%），约 0.5 h |
| 62 | 1e-3 | stretch → 512, Nv = 1000 → 2000 | 0.0005 | 342800 | 同 Case59，Nv = 2000：检验 4 τi 时壁面 Te 是否受 Nv = 1000 速度网格限制，约 1.8 h（均匀 4000 × 2000 约 13 h） |

Case58–60 与 Case50 只差 x 网格（壁面格宽相同 0.5λ），用于检验"只在壁面加密"能否复现标杆结果、体区最粗能用多少。
Case61、62 应在 Case59 与 Case50 一致之后再看（若体区 8λ 不够，按 Case58 的网格重做）。

Case59 结果：0–4 τi 壁面粒子流 / 能流与 Case50 相差 ≤ 0.11%，T_i·t² 相同，壁面 Te 相差 ≤ 0.5%；工作站用时 58 min（Case50 约 6.7 h）。

| Case | λ | Nx | Δt | 步数 | 说明 |
|---|---|---|---|---|---|
| 63 | 1e-3 | refine → 280, Nv = 1000 | 0.0005 | 342800 | `XGrid refine`：均匀粗网格 8λ（250 格），每侧壁面 1 个粗格等分为 16 个 0.5λ 细格（细网格 8λ），格宽在粗细界面处突变 16 倍，约 32 min |
| 64 | 1e-3 | refine → 310, Nv = 1000 | 0.0005 | 342800 | 同 Case63，`NRefine 2`（细网格 16λ），约 35 min |
| 65 | 1e-3 | refine → 370, Nv = 1000 | 0.0005 | 342800 | 同 Case63，`NRefine 4`（细网格 32λ），约 42 min |
| 66 | 1e-3 | refine → 430, Nv = 1000 | 0.0005 | 342800 | 同 Case63，`NRefine 6`（细网格 48λ）：与 Case59（50λ 细网格 + 1.05 渐变）只差有无过渡区，约 49 min |

Case63–66 用于检验"粗网格 + 壁面细网格块"（便于移植 kobra）：细网格区需要多宽，以及格宽突变是否引入误差。

Case63 结果：每侧只细分 1 个粗格（细网格 8λ，16 倍格宽突变）即与 Case50 吻合：0–4 τi 壁面粒子流 / 能流相差 ≤ 0.08%，T_i·t² 相同，
壁面 Te 相差 ≤ 0.5%，与 Case59 几乎相同；工作站用时 21 min（Case50 约 6.7 h，Case59 58 min）。

以 Case63 为基准（refine 网格：粗格 8λ、细格 0.5λ）的扫描，Case61、62 被 Case67、72 取代，不必运行：

| Case | λ | Nx | Δt | 步数 | 说明 |
|---|---|---|---|---|---|
| 67 | 1e-3 | refine → 280, Nv = 1000 | 0.001 | 171400 | 同 Case63，Δt = 0.001（细格 x-CFL = 12），约 11 min |
| 68 | 1e-3 | refine → 280, Nv = 1000 | 0.002 | 85700 | 同 Case63，Δt = 0.002（细格 x-CFL = 24，ω_pe Δt = 2），约 6 min |
| 69 | 1e-3 | refine → 264, Nv = 1000 | 0.0005 | 342800 | 同 Case63，`DxWall 1`：壁面粗格等分为 8 个 1λ 细格（细格分辨率下限：0.5λ 已分辨鞘层，2λ 的均匀网格 Case48/49 后期升温），约 20 min |
| 70 | 1e-3 | refine → 256, Nv = 1000 | 0.0005 | 342800 | 同 Case63，`DxWall 2`：壁面粗格等分为 4 个 2λ 细格；与 Case49（均匀 2λ）只差体区网格，预计重现后期升温，约 19 min |
| 71 | 1e-3 | refine → 280, Nv = 500 | 0.0005 | 342800 | 同 Case63，Nv = 500：速度分辨率下限，约 11 min |
| 72 | 1e-3 | refine → 280, Nv = 2000 | 0.0005 | 342800 | 同 Case63，Nv = 2000，约 42 min |
| 73 | 1e-3 | refine → 280, Nv = 4000 | 0.0005 | 342800 | 同 Case63，Nv = 4000：Nv 收敛参照（Case50 的壁面 Te 在 t ≳ 2.5 τi 可能受 Nv = 1000 限制），约 1.4 h |
| 74 | 1e-3 | refine → 280, Nv = 1000 | 0.0005 | 342800 | 同 Case63，`SolverType Poisson`（标准泊松，ω_pe Δt = 0.5），与 Case54（均匀 4000 格）对照，约 21 min |
| 75 | 1e-3 | refine → 280, Nv = 1000 | 0.001 | 171400 | 同 Case67，`SolverType Poisson`（ω_pe Δt = 1）：判断 Δt = 0.001 时壁面 Te 偏低来自 AP 项还是时间离散，约 11 min |
| 76 | 1e-3 | 250（均匀 8λ）, Nv = 1000 | 0.0005 | 342800 | 同 Case74，但均匀 8λ 网格（无壁面细格）：看普通泊松在粗网格上失稳的时间、位置是否与 Case74 相同；OutF 8570（每 0.1 τi 输出一次 f），约 19 min |
| 77 | 1e-3 | refine → 280, Nv = 1000 | 0.0005 | 51420 | 同 Case74，只算到 0.6 τi，OutF 1714（每 0.02 τi 输出一次 f）：确定普通泊松的失稳从细格内部还是粗细界面（8λ）开始，约 4 min |
| 78 | 1e-3 | refine → 280, Nv = 1000 | 0.0005 | 342800 | 同 Case74（refine + 普通泊松），`WallE onesided`：普通泊松下矩阵不变，只改贴壁格 E₀ = −(φ₁ − φ₀)/Δx，看贴壁格电场是否为失稳的触发点；`OutProf 1714`（剖面每 0.02 τi），约 21 min |
| 79 | 1e-3 | 250（均匀 8λ）, Nv = 1000 | 0.0005 | 342800 | 同 Case76（均匀 8λ + 普通泊松），`WallE onesided`，`OutProf 1714`，约 19 min |
| 80 | 1e-3 | refine → 280, Nv = 1000 | 0.0005 | 342800 | 同 Case63（refine + AP），`WallE onesided`：Case78 的 AP 对照（onesided 本身对结果的影响），`OutProf 1714`，约 21 min |
| 81 | 1e-3 | stretch → 512, Nv = 1000 | 0.0005 | 342800 | 同 Case59（壁面 50λ 内 0.5λ，等比 1.05 过渡到 8λ），`SolverType Poisson` + `WallE onesided`：粗细之间无突变时普通泊松是否稳定（对照 Case78 的 16 倍突变），`OutProf 1714`，约 1 h |
| 82 | 1e-3 | 250（均匀 8λ）, Nv = 1000 | 0.0005 | 342800 | 同 Case79，`SolverType AP`：全粗网格上 AP 与普通泊松的精度对照，约 19 min |

Case67–69 结果（以 Case50 为基准）：
- Case67（Δt = 0.001）与均匀网格 Case52 逐时刻几乎相同（壁面 Te 低 2–5%，通量 < 0.1%）。
- Case68（Δt = 0.002）稳定，通量与 T_i·t² 仍 < 0.1%，壁面 Te 偏低更多：2 τi −10%，4 τi −14%。
- Case69（细格 1λ）失败：壁面 Te 在 1 τi 时 +36%，1.5 τi 时 +100%（数值升温），细格须 ≤ 0.5λ。
- Case74、75（refine 网格 + `SolverType Poisson`，Δt = 0.0005 / 0.001）失败：0.4 τi 时与 Case63 相差 0.1%，约 0.45–0.5 τi 起壁面约 20λ 内失稳
  （0.8 τi 贴壁格电子密度 −93%，电场能约 100，正常约 0.1），体区不受影响，减小 Δt 无效。均匀 0.5λ 网格上的普通泊松（Case54）正常，
  故原因为体区粗格 Δx ≫ λ_D,loc（0.5 τi 时约 25 倍）：普通泊松要求 Δx ≲ λ_D，AP 格式无此要求。非均匀 / AMR 网格必须用 AP。
- Case76（均匀 8λ + 普通泊松）：同样在壁面失稳（贴壁格电子被抽空、格点尺度锯齿），但更早、更强：0.2 τi 出现锯齿，0.3 τi 电场能偏离，
  饱和约 470（Case74 约 90），体区仍不受影响。即失稳由鞘层电荷分离区与粗格接触触发；壁面细网格块只推迟、减弱，不能消除。
- Case77（Case74 每 0.02 τi 输出 f，与 Case74 逐位相同）：电场能 0.46 τi 开始偏离，0.48→0.52 τi 爆发。0.50 τi 时电子平均速度的扰动覆盖离壁 0–40λ，
  在 8–16λ 附近变号，16–40λ（紧接细网格块的前几个粗格）幅度最大；> 40λ 的体区不受影响。密度变化最先在贴壁格（0.25λ）显现（该处密度最低），
  随后贴壁格电子在 v > 0 / v < 0 之间来回振荡（壁面电场反复翻转）。即失稳模式跨过粗细界面，位于鞘层 / 预鞘层与粗格相接处，而非体区。
- Case78（Case74 + `WallE onesided`）同样失败，时间与 Case74 相同：电场能 0.50 τi 偏离（×1.9），0.52 τi ×27。0.46–0.48 τi 时 8–40λ（细网格块之后的粗格）
  电场先变大并反号（|E| 0.1 → 0.6 → 2.7），贴壁 0–8λ 的鞘层电场（约 11）几乎不变；失稳前 Δx/λ_D,loc 在细格约 1、粗格 12–15。
  即贴壁格电场的离散方式不是触发点，失稳起于预鞘层中 Δx ≫ λ_D,loc 的粗格。失稳前电场能比 Case74 低约 15%，是因为 onesided 去掉了贴壁格从格中心到壁面的电势降
- Case79（均匀 8λ + 普通泊松 + `WallE onesided`，到 1.65 τi）稳定：电场能平滑下降、无失稳（同一网格 central 的 Case76 在 0.3 τi 失稳）。
  与 Case63 比，壁面粒子流 / 离子能流 < 0.05%，电子能流 < 0.5%，壁面 Te 0.5 / 1 / 1.5 τi 时 +5.5% / +3.6% / −2.2%。
  电场能只有 Case63 的 2–20%（8λ 贴壁格内鞘层不分辨，onesided 下鞘层电场不作用在贴壁格上）。
  对比 Case78（refine + onesided 失败）：失稳与粗细格 16 倍突变处的普通泊松离散有关，而非粗格本身；Case81 检验平滑过渡时是否稳定

Nv 的经验判据（Case46 Nv = 256 对 Case50 Nv = 1000，同 Nx = 4000、Δt）：σ/Δv（σ = √T_wall 为壁面速度宽度，Δv = 2·Umax/Nv）
降到约 1.5 时离子 T 误差约 1%、降到约 0.75 时约 5%；电子 Te 对此更敏感，σ/Δv ≈ 3.6 / 1.7 / 1.05 时误差约 1% / 4% / 10%。
壁面粒子流、能流对 Nv 不敏感（Nv = 256 到 2.75 τi 误差仍 < 0.3%）。

Case52–54 结果（Nx = 4000，Nv = 1000，以 Case50 为基准；Case53、54 为中途数据，分别到 2.85 / 3.05 τi）：
壁面粒子流、能流（离子与电子）与 Case50 相差 < 0.01%，T_i·t² 相同（4 τi 时 0.00505），Γe/Γi − 1 相同；
只有壁面 Te 随 Δt 单调变化：Δt = 0.001 低 5.2%（4 τi），Δt = 0.00025 高 2.2%（2.5 τi）。标准泊松（Case54）与 AP 相差 ≤ 0.6%。
即鞘层分辨后，Δt 与 AP 项对结果基本无影响。

Case51–53 的 OutDiag / OutF 按 Δt 调整，使诊断间隔（约 0.004–0.005）和分布函数输出次数（10 次）与 Case50 相同。
Δt = 0.004（ω_pe Δt = 4，x-CFL = 48）在 Nx = 4000 下于 t ≈ 0.008 τi 发散（见下），不列入扫描。

Case46–48 的 Δt 统一取 0.0005（x-CFL = 1.5 / 6 / 12）：Nx = 4000 时 Δt = 0.004（x-CFL = 48）在电子首次到达壁面（t ≈ 0.008 τi）时电场发散、粒子全部丢失，
与幽灵格层数无关（nghosts = 6 结果逐位相同）；x-CFL ≤ 24 时可通过该时刻。旧 Case14 的结果即为 Δt = 0.004 下发散的数据，作废。

## case 文件可选参数（缺省值）
- `SolverType`（必填）：`AP`（渐近保持重构泊松）或 `Poisson`（标准泊松 λ²φ'' = -ρ：同一 Thomas 求解器、壁面离散与电场计算，只去掉 Δt²n/m、Δt∂J、Δt²∂²R 项；显式格式要求 ω_pe Δt < 2 且 Δx ≲ λ）。旧的 `SOR` / `SORAP` / `DS` / `Thomas` 未维护，不要使用
- `lambda`（1e-3）：λ_D / L_par
- `XGrid`（uniform）：`uniform` 为 Nx 等分；`stretch` 为非均匀 x 网格（壁面细、体区粗，Nx 由下列参数算出，case 中的 Nx 不用）。x 输运为非均匀网格上的守恒半拉格朗日 PFC3（全局 Δt，可跨越多格），AP / 泊松方程为非均匀有限体积；只支持 Scheme PFC3、SolverType AP / Poisson、WallE central / onesided、WallMoment avg、EFilter 0
- `DxWall`（0.5）、`WallWidth`（50）、`DxBulk`（16）：stretch 网格的壁面格宽、壁面等宽区宽度、体区格宽，单位 λ
- `Stretch`（1.05）：stretch 网格过渡区相邻格宽之比
- `XGrid refine`：均匀粗网格（格宽 `DxBulk` λ，要求 2/(DxBulk·λ) 为整数），每侧壁面最近的 `NRefine`（1）个粗格各等分为 DxBulk/DxWall（须为整数）个细格，细格面与粗格面重合、无过渡区；Nx = Nc + 2·NRefine·(r − 1)。`NRefine 0` 与均匀网格逐位相同。数值方法与 stretch 相同
- `WallFace`（PFC3）：`PFC3` 或 `upwind1`（要求 x-CFL < 1）
- `VCFL`（0）：v 方向子循环的 CFL 上限，0 为不子循环
- `OutDiag`（10）：`runtime_diags.csv` 输出间隔（步）
- `OutF`（100）：分布函数快照输出间隔（步），0 为不输出
- `OutProf`（OutF/10）：剖面 `out/prof.vpd` 输出间隔（步），0 为不输出
- `OutPrec`（8）：vpd 快照中 f 的字节数，8（double）或 4（float，文件减半）
- `EFilter`（0）：内部电场 (1,2,1)/4 二项式滤波次数，贴壁 2 格不滤
- `SpikeDump`（0）：电场能突跳诊断。某步电场能超过上一步的 `SpikeDump` 倍时，在 `spike/ev<k>_` 下保存触发步的 fⁿ 与前一步 f^{n-1}（`fe_now/fe_prev/fi_now/fi_prev.txt`，文本，nv 行 × nx 列），以及触发前后各 20 步的场剖面 `prof.csv`（φ, E, ρ, J, R, Σf/m，含 i = −1 与 Nx 两个幽灵点）；最多 5 次
- `WallE`（central）：贴壁格电场离散。`central`：原方案，E₀ = −(φ₁ − φ_g)/2dx，φ_g = −φ₀；`onesided`：E₀ = −(φ₁ − φ₀)/dx，且 AP 方程壁面面系数只保留 λ²（满足离散高斯定律，后期电子数值加热大幅减小）；`kobra`：与 quasi-kobra 的 mug_kobra.c 相同（幽灵 φ_g = −2.5φ₀ + 0.5φ₁，壁面电流取壁面面通量，壁面 ∂R 单侧估计）
- `WallMoment`（avg）：只对 `WallE central` 起作用，AP 方程壁面面上的矩。`avg`：原方案，壁面面系数 λ² + Δt²·½(n_幽灵 + n₀)、面电流 ½(J_幽灵 + J₀)、∂R 用幽灵格中心差分（幽灵格只含出射方向的 f，相当于 n_out + ½n_in、J_out + ½J_in）；
  `out`：与 kobra 相同，面系数 λ² + Δt²·n_out/m、面电流 J_out、∂R 取入射方向单侧估计；φ 幽灵格（φ_g = −φ₀）与贴壁格 E₀ 不变
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

画图（新算例 `out/*.vpd`，旧算例的 txt 也能读；离壁距离以 λ 为单位，缺省两壁平均、矢量取指向壁面为正，虚线为粗细界面）：
- `~/plotenv/bin/python plot/VPplot.py prof runs/Case78 --t 0.3,0.46,0.5 --vars E,phi,dn,n_e,T_e,dxlam --dmax 200 --logd`：剖面（E、φ、n_i−n_e、n、u、T、q、局部 λ_D、Δx/λ_D,loc），多时刻用由浅到深的蓝色；
  给多个算例时颜色区分算例、线型区分时刻
- `... VPplot.py xt runs/Case78 --var E --dmax 200 --logd --tmin 0.3 --tmax 0.6 --vmax 3`：离壁距离 × 时间的时空图；`--ref runs/Case80` 画与对照算例之差
- `... VPplot.py f runs/Case78 --t 0.6 --sp e --dmax 60 --log`：壁面附近相空间 f(x, v)（v 范围自动，`--vlim` 指定）
- `... VPplot.py fv runs/Case78 runs/Case80 --t 0.6 --d 0.25,10,100 --log`：给定离壁距离处的 f(v)
- `... VPplot.py diag runs/Case63 runs/Case78 --cols ElectricEnergy,Tw.e- --logy`（`--rel` 为相对第一个算例）
- 分布函数（`plot/VPdist.py`，离子速度以 v_Ti 为单位，色标缺省对数）：
  - `... VPdist.py panel runs/Case79 --t 0.2,0.5,1.0,1.6 --dmax 200 --logd`：相空间小多图（行 = 电子 / 离子，列 = 时刻，同一粒子共用色标）
  - `... VPdist.py anim runs/Case79 --dmax 200 --logd`：全部快照的相空间动画（GIF）
  - `... VPdist.py slice runs/Case79 runs/Case63 --t 0.4,1.2 --d 0.5,20,200 --sp e --maxw`：给定离壁距离的 f(v)，`--maxw` 叠加同 n、u、T 的麦克斯韦分布
  - `... VPdist.py vt runs/Case78 --d 12 --sp e`：给定离壁距离处 f(v) 随时间的演化（时间分辨率为 OutF）
  - `... VPdist.py diff runs/Case78 runs/Case80 --t 0.6 --dmax 100`：两个算例同一时刻 f 之差（同一网格）
- 图默认存在第一个算例目录下，`--out` 指定文件名，`--side left|right|full` 只看一侧或全域；`python plot/VPplot.py -h` 看全部选项
- `~/plotenv/bin/python plot/PaperFigs.py runs/Case4 [runs/Case12 ...]`：按文献坐标轴画图 1、图 7 和 Γe/Γi − 1，图存在运行目录（多个目录时存在 runs/）

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
