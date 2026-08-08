# cpp_lite — 冻结分支后的精简版 Fourier_Quad pipeline (C++)

本目录是 `../cppv2` 的精简重构版，与 `../f77_lite` 一一对应：把 8 个编译期分支
开关固定为当前使用的取值，并**物理删除**所有未选中分支的代码和相应的 `if` 语句。

当前集成结构中，数值实现统一位于 `src/process_main/`，对应头文件位于
`include/process_main/`；初始化器独立位于 `src/process_init/` 与
`include/process_init/`。下文未写目录的数值文件名均指这两个 `process_main` 目录。

未修改的文件（`FitsIO.*`、`ImageProcessing.*`、`LinearSolve.*`、`NumericalRecipes.*`、
`UniversalUtils.*`、`ExStar.*`、`ExposureInfo.*`、`FourierTransformSt2.*`、
`FourierTransformSt1.hpp`、`CatalogCombiner.hpp`、`ShearMeasurement.hpp`）
是逐字节拷贝，未作任何改动。

## 一、冻结的开关

| 开关 | 冻结值 | 保留下来的行为 |
|---|---|---|
| `ASTROMETRY_trivial` | 0 | 只走 Gaia 星表天测（`genAstrometryData` / `getAstrometry`） |
| `include_FLAT` | 0 | 不做 super-flat 乘算 |
| `include_Mask` | 2 | 逐 chip 读 `dirOutput/dqmask/<expo>_<cid>.fits` 掩膜 |
| `ext_cat` | 1 | 使用外部星表 `SOURCE_CAT` 建源目录 |
| `ext_PSF` | 0 | PSF 由本帧恒星测量 |
| `deblending` | 1 | 恒定调用 `deBlending` |
| `PSF_type` | 1 | 局域多项式 PSF 拟合（`makePSFLocalFit`） |
| `PSF_Ms` | 0 | 不做多尺度 / PCA PSF 重建 |

## 二、仍然保留的可选分支

`PROCESS_stage`、`CCD_split`、`gal_smooth`、`star_smooth` 的分支**原样保留**，
常量仍在 `include/process_main/LensingConfig.hpp` 中。

## 三、逐文件改动

| 文件 | 行数 | 改动 |
|---|---|---|
| `LensingConfig.hpp` | 177 → 170 | 删 8 个开关，以及 `FLAT_PATH`、`PSF_PATH`、`flat_thresh`、`step_psf`、`n_neighbor` 与整个 PCA 参数段（`rescale_size` `procs_pn` `work_pn` `nblocks` `n_pcs` `npp6th` `pca_negative_eigenvalue_threshold` `nmax_star_pchip`）；F6 mode-bar 参数与主树保持一致 |
| `MPIScheduler.hpp` | 22 → 21 | 删 `forcecov()` 声明（只被 PSF_Ms 的 PCA 阶段使用） |
| `MPIScheduler.cpp` | 163 -> 91 | 新增（原 `cppv2/` 缺失，已补回）；删 `forcecov()` 定义，与已重构的 `.hpp` 对齐 |
| `Makefile` | 43 → 42 | 从 `SRCS` 移除 `PSFRecons.cpp` |
| `main.cpp` | 203 → 195 | 删 `PSFRecons.hpp` include、`chipPSFRecons`、`freePSFMemory` 调用 |
| `PreProcess.hpp/.cpp` | 1526 → 1471 | `chipPreProcess` 去掉 `maskFile` 形参；删 flat 读取与乘算；掩膜段直接采用 `include_Mask=2` 实现；天测段直接走 Gaia 分支；当前 F6 核心与主树一致 |
| `Astrometry.hpp/.cpp` | 1687 → 1617 | 删 `genAstrometryDataTrivial`、`getAstrometryTrivial`；`chipProcessAstrometry` 直接调 `getAstrometry`；恢复本地自包含头文件 |
| `SourceExtractor.hpp/.cpp` | 1253 → 908 | `chipProcessSource` 去掉 `flatFile` 形参；删 `getFlatName`、`genSourceCatalog`、`genStarCandidate`；`ext_cat=1` 分支展平；增加失败 norm chip 入口守卫 |
| `PSFModel.hpp/.cpp` | 1544 → 942 | 删 `initAndLoadAllPSF`、`freePSFMemory`、全局 PCA 存储与索引 helper、`PSF_rescale`、`PSF_unscale`、`saveRescaleFactor`、`makePSFHybrid`、`interpolatePSF`、`getPSFModelVeryLocal`、`genPSFFits`；`makePSFLocalFit` 内所有 `PSF_Ms=1` 段删除 |
| `PSFRecons.cpp/.hpp` | 975 → 删除 | 整个文件只服务于 PSF_Ms=1 |
| `ShearMeasurement.cpp` | 490 → 425 | 删 `PSFRecons.hpp` include、外部 PSF 读取段、rescale factor 读取段、`PSF_type=2` 段；PSF 模型求值直接走 `PSF_type=1 / PSF_Ms=0` 路径 |
| `FourierTransformSt1.cpp` | 107 → 103 | 删 `if (ext_PSF == 1) return` |
| `CatalogCombiner.cpp` | 216 → 166 | 两处 `ext_cat == 1` 分支展平，删 `else` 分支 |

当前 `process_main` 数值 `.cpp/.hpp` 合计 **12347 → 10071 行**；新增 F6 mode-bar 核心在主树与 lite
中逐字一致，不改变冻结分支的行为差异。

另外，`PSFModel.cpp` 内部静态辅助函数 `fitPSFCoefficients` 的
`bool normalize_positions` 形参与对应的 `if (normalize_positions)` 一并删除：
该开关唯一的 `false` 传入方是 `interpolatePSF`（PSF_type=2 的 hybrid 拟合），
删除后仅剩 `itpNormPSF` 恒传 `true`，故坐标归一化改为无条件执行。

## 四、刻意保留的东西

`LinearSolve.cpp/.hpp`、`UniversalUtils.cpp/.hpp`、`FitsIO.cpp/.hpp`
**未做任何改动**。它们是通用数值 / IO 工具库，地位等同 f77 树里的
`FFTPACK.f` / `press.f`，因此即使个别函数随 `PSFRecons` 一起失去了调用者，
也整体保留。失去调用者的清单（供你决定是否进一步清理）：

| 函数 | 原调用者 |
|---|---|
| `LinearSolve::analyzeCovarianceSpectrum`（及 `EigenSpectrumDiagnostics`） | `PSFRecons::chipResPCAFit` |
| `UniversalUtils::fit2D2` | `PSFRecons::itpNormPSFCov` |
| `UniversalUtils::fit2D2Cov`（它本身是 `LeastSquaresQR::unscaledCovariance` 的唯一调用者） | `PSFRecons::itpNormPSFCov` |
| `UniversalUtils::invertMatrixF77` | `PSFRecons::fitPcaBlockCoefficients` |
| `FitsIO::readPara` | `PSFModel::makePSFHybrid`（f77 侧 `read_para` 同样失去调用者） |

注意 `LeastSquaresQR::unscaledCovariance` **仍被调用**（来自
`UniversalUtils::fit2D2Cov`），它不是孤儿。

以下是**原版就已经存在**的死代码，与本次重构无关，未作处理：
`ShearMeasurement::getWindowMinKVer2`（f77 侧 `get_window_min_k_ver2` 同样从未被调用）。

`PSFModel.cpp` 中 `makePSFLocalFit` 里的 `int nplx = LensingConfig::nplx;`
在原版就是未使用变量（`itpNormPSF` 的签名不含 `nplx`），本次一并保留，
`-Wall` 下的警告与原版相同。

`PreProcess.cpp` 的掩膜错误信息仍沿用原版措辞
（`"Error / wrong size of flat file!"`、局部变量名 `flat_weight`）。这段字符串
在原版的 `include_Mask==2` 分支里就是这样写的，改动会让程序输出与原版不一致，
故保持原样——但请注意这里指的是 DQ mask，不是 super-flat。

`CatalogCombiner::combineExpoCatalog` 里 `g1c`/`g2c` 在保留分支中恒为 0，
两行 shear 修正因此是恒等运算；`g1_c`/`g2_c` 常量只剩注释引用。这与原版
`ext_cat=1` 路径完全一致（原版就是这样写的），保留以便随时恢复该修正。

失去引用的常量（保留在 `include/process_main/LensingConfig.hpp`，供你决定是否清理）：
`Camera_ccd_num`、`psf_order`（原版就未使用）、`npox`（原版就未使用）、
`g1_c`、`g2_c`（仍被 `CatalogCombiner.cpp` 中注释掉的代码引用）、`nplx`。

不再被写入的输出目录（也确认没有任何 lite 代码去**读**它们）：
`rescale/`、`starxy/`、`fits_psfresi/`、`dat_pcs/`、`dat_starcomp/`。

## 五、原本缺失的 `MPIScheduler.cpp`（已补回）

`Makefile` 的 `SRCS` 里列了 `MPIScheduler.cpp`，但该文件此前不在 `cppv2/` 目录中，
`MPIScheduler::init/finalize/barrier/distribute/my_id/num_procs` 全部无定义，原始树
同样无法链接。现已补回：`src/process_main/MPIScheduler.cpp` 按重构规则删去
`forcecov()` 定义（只被已删除的 `PSFRecons` 使用，与已重构的 `MPIScheduler.hpp`
对齐），其余 `init/finalize/barrier/distribute` 原样保留。上述符号全部解析正常。

## 六、编译验证

```bash
make CXX=mpicxx STACK_PREFIX=<scientific-stack-prefix> \
  EIGEN_INCLUDE=<eigen3-include-directory> -j4
```

16 个数值 `.cpp` 源文件以及统一入口、`process_main` wrapper 和三个
`process_init` 目录实现
全部通过编译，零错误；与原始 `cppv2` 数值树逐文件对照，精简树仅少了被删的
`PSFRecons.cpp`。`-Wall` 警告与原版相同且无新增（`Astrometry.cpp`
2 条、`PSFModel.cpp` 2 条，均为原版就有的 unused-variable；`MPIScheduler.cpp` 零警告）。
完整链接成功，此前未解析的 `MPIScheduler::*` 已由补回的 `MPIScheduler.cpp` 提供。
`.o` 层未定义符号集合（扣除 cfitsio / FFTW / LAPACK / BLAS / MPI 库符号）为原始树的
真子集：lite 158 个、原始 164 个，**零新增**未解析符号。

2026-07-30 再验证：本地 `Astrometry.hpp` 已恢复，lite 不再需要 `-I../cppv2`；完整
C++17 编译和链接成功。新增 `tests/SetSigVerification.cpp` 输出 `PASS`，并与主 C++
得到相同结果：纯噪声恢复 `1.001715`、拥挤度漂移 `0.2580%`、归一化 RMS
`1.002392`，失败 amp 不修改自身图像且不回滚已应用的前一 amp。
