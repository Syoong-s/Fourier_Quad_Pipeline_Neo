# Fourier_Quad C++ 流水线指南

本文说明如何编译、配置和运行 `Fourier_Quad_Cpp`。

> English: [CPP_GUIDE.md](CPP_GUIDE.md)

## 程序结构

[`cpp_Standard`](cpp_Standard/) 与 [`cpp_Lite`](cpp_Lite/) 均生成
`Fourier_Quad_Pipe`，并包含：

- `main.cpp`：MPI 生命周期与六个顶层阶段的调度；
- `pipeline.example.ini`：完整运行期配置模板；
- `config/`：编译默认值与固定数值参数；
- `include/`、`src/`：阶段模块、运行配置与星表布局；
- `include/general/`、`src/general/`：曝光表、路径、MPI、调度与通用数值工具；
- `Makefile`：C++17/MPI 构建入口。

Standard 保留平场、掩膜、简化测天、外部/混合/PCA PSF 等可选分支。Lite 只保留
Gaia 测天、逐 CCD DQ 掩膜、外部源星表、帧内恒星、去混叠、局域多项式 PSF，且
不含 PCA。Lite 是物理删除未使用分支后的版本，不是换了默认值的 Standard。

## 处理流程

六个顶层阶段按固定顺序执行：

| 阶段 | CLI | 作用 |
|---|---|---|
| `process_astrocat` | `--run-astrocat` | 将原始两列 Gaia 星表重分块、去重为一度瓦片。 |
| `process_extcat` | `--run-extcat` | 将原始 External source catalog 文件重分块为天空瓦片。 |
| `process_init` | `--run-init` | 搜索 `.fits.fz`，提取 Science/DQ CCD，生成曝光表。 |
| `process_main` | `--run-main` | 执行九阶段数值流水线。 |
| `process_rearr` | `--run-rearr` | 将 `*_all.cat` 空间分区。 |
| `process_fd` | `--run-fd` | 按场畸变分箱恢复平均剪切。 |

`process_astrocat`、`process_extcat` 每次运行依次各执行一次，其余阶段按数据集顺序执行；
首次集体失败会停止整个任务。

`process_main` 用素数乘积选择阶段：

| 阶段 | 素数 | 内容 |
|---:|---:|---|
| 1 | 2 | 背景/噪声预处理与 Gaia 匹配 |
| 2 | 3 | 测天解 |
| 3 | 5 | 源检测、去混叠、恒星候选体 |
| 4 | 7 | 恒星候选体功率谱 |
| 5 | 11 | PSF 选择与建模 |
| 6 | 13 | 星系功率谱 |
| 7 | 17 | Fourier_Quad 估计量与形态测量 |
| 8 | 19 | 曝光级统计 |
| 9 | 23 | 星表合并与标定 |

`process_stage` 能被某素数整除时执行对应阶段。完整值 `223092870` 启用全部阶段；
阶段 9 必须与阶段 8 同时启用。

## 编译

需要支持 C++17 的 MPI C++ 编译器，以及 CFITSIO、FFTW3、Eigen3、LAPACK、BLAS。

```bash
cd cpp_Lite                    # 或 cpp_Standard
make -j4
./Fourier_Quad_Pipe --help
```

科学库集中安装在同一前缀时：

```bash
make CXX=/path/to/mpicxx STACK_PREFIX=/opt/science-stack -j4
```

Eigen 单独安装时再传入 `EIGEN_INCLUDE=/opt/eigen/include/eigen3`。当前 Standard 与
Lite Makefile 只公开 `all` 和 `clean`；`make` 等同于 `make all`。修改编译期配置或
工具链后，先执行 `make clean` 再重新编译。

Stage 5 使用所有同 CCD FWHM-locus 配对计算每颗候选体的 minChi，但曝光阈值只由
触及受限大尺寸 reference 的唯一配对估计。raw analytic PRESS 保留为诊断；可选 rejection
使用 leverage-standardized PRESS 和受保护的临时重拟合。关闭 rejection、触发删除保护或
重拟合失败时均保留合法 first fit。这些科学开关仍是编译期设置，已有运行时 PSF 模式、
芯片几何与 direct stamp-cube I/O 不变。

## 配置

先复制所选版本的模板：

```bash
cp pipeline.example.ini pipeline.ini
```

优先级为：

```text
编译默认值 < INI < CLI
```

INI 分为五段：`[process]` 控制阶段和输出路径，`[astrocat]` 控制原始两列 Gaia 输入、
独立输出目录、header 与已有输出策略，`[extcat]` 控制原始星表解析与布局，`[init]`
控制归档根目录和数据集，`[lensing]` 控制可在运行间改变的科学分支与相机参数。
Lite 会拒绝 Standard 专属键；未知段、未知键、无效值及不一致的阶段/schema 组合均会在
执行前报错。

`[lensing].astrometry_cat_type=1` 读取旧式大 Gaia 瓦片；值 `2` 累积读取
`process_astrocat` 生成的一度瓦片。生产者路径 `[astrocat].output_directory` 与消费者
路径 `[lensing].astrometry_cat` 始终是两个独立设置；若同一次运行的后续阶段需要消费
新发布的瓦片，应显式分别配置两者。

CLI 同时接受 `--name value` 与 `--name=value`。布尔值接受 `true/false`、`1/0`、
`yes/no`、`on/off`。首次显式 `--dataset`、`--contains`、`--extcat-contains` 会替换
INI 列表，后续同名选项追加。一个裸位置参数仍可作为 `--expo-list` 的兼容写法。

完整 CLI 以 `./Fourier_Quad_Pipe --help` 为准；所有 INI 键见
[CPP_PIPELINE_PARAMETERS.md](CPP_PIPELINE_PARAMETERS.md)。

运行期设置优先写入所选版本的 INI，CLI 只用于单次覆盖。头文件专属设置应修改
`cpp_Standard/config/` 或 `cpp_Lite/config/` 下对应的配置头文件；修改源参数而非
派生常量，然后执行 `make clean` 并重新编译。参数参考按七个配置头文件各给一张
完整的 Standard/Lite 对照表，并明确标出哪些值可在不重编译的情况下覆盖。

## 常用运行方式

按 INI 串联运行：

```bash
mpirun -np 8 ./Fourier_Quad_Pipe --config pipeline.ini
```

仅对已有曝光表运行主流程：

```bash
mpirun -np 8 ./Fourier_Quad_Pipe --config pipeline.ini \
  --run-init false --run-main true --run-rearr false --run-fd false \
  --expo-list /data/work/expo_gband.list
```

仅重分块原始 Gaia 星表：

```bash
mpirun -np 8 ./Fourier_Quad_Pipe \
  --run-astrocat true --run-extcat false --run-init false --run-main false \
  --run-rearr false --run-fd false \
  --astrocat-input /data/raw_gaia --astrocat-output /data/gaia/tiles \
  --astrocat-add-header true --astrocat-existing fail
```

`--astrocat-output` 只控制 `process_astrocat`，不会与
`[lensing].astrometry_cat` 比较、同步，也不要求两者一致。

仅生成外部星表瓦片：

```bash
mpirun -np 8 ./Fourier_Quad_Pipe \
  --run-extcat true --run-init false --run-main false \
  --run-rearr false --run-fd false \
  --extcat-input /data/raw_catalogs --extcat-output /data/catalogs/tiles
```

至少启用一个顶层阶段。初始化成功后，后续阶段自动使用生成的绝对路径
`expo_<target>.list`。在 Slurm 上应改用站点支持的启动器；本仓库容器 runner 固定使用
`srun --mpi=pmi2`。

## 输入与输出

Science images、Gaia catalog、External source catalog 与 DQ masks 的统一定义和最低
要求见根目录 README 的[输入数据要求](README_CN.md#输入数据要求)。DQ masks 是取决于
配置的输入类别：Lite 始终读取逐 CCD DQ masks；Standard 只有在所选
`[lensing].include_mask` 模式和实际执行路径都不访问 DQ 数据时才可省略。

曝光表每个非空记录包含一个 CCD 列表路径，可带兼容用的 CCD 数量；当前 C++ 读取器
支持带引号的路径。

初始化器读取原始压缩归档但不移动它们，在 `output_root/<target>/` 下创建 `science/`、
`dqmask/`、`stamps/`、`result/`，并发布 `expo_<target>.list`、
`fits_<target>.list` 与 manifest。主要下游结果为：

```text
<dataset>/result/<exposure>_all.cat
<dataset>/<rearr-output-dir>/subcat_*.cat
<dataset>/<rearr-output-dir>/catalog_summary.txt
<dataset>/<fd-output-dir>/FD_test_comb.dat
```

阶段 7 输出 28 个流水线字段；阶段 9 写入 `EXPO_NUM`、`ccD_NUM`，并附加一个曝光
`chi2`。默认完整行因此是 18 个外部字段、两个身份字段和 29 个流水线字段，共 49 列。
显式投影只改变外部前缀宽度；身份顺序始终是 `EXPO_NUM` 紧邻并位于 `ccD_NUM` 之前，
FD 读取该真实曝光身份，不再从文件列表顺序推断。旧的 48 列产物与新 schema 不兼容，
在重排或 FD 前必须重新生成。列身份、可选星等映射与 FD 波段选择顺序见
[CPP_PIPELINE_PARAMETERS.md](CPP_PIPELINE_PARAMETERS.md)。

无效数值源仍保持行号对齐：阶段 6 写 12 个 `-99999` 标记，阶段 7 对该源或后续非有限
结果写完整 28 个 `-99999`，阶段 9 会丢弃该行。这是输出约定，不是额外配置项。

## 常见错误

- 不要在未启用阶段 8 时启用阶段 9。
- 外部星表输出目录不能等于或位于输入目录之下。
- 显式投影必须包含原始 RA、Dec、photo-z 列。
- Lite 不能使用 Standard 专属 `[lensing]` 键。
- 容器内程序参数必须使用容器路径，而不是宿主路径。
- 不要让多个作业同时清理或编译同一源码副本。

本地 Docker 见 [cpp_docker/README-CN.md](cpp_docker/README-CN.md)，Slurm/Apptainer
见 [cpp_docker/runner/README-CN.md](cpp_docker/runner/README-CN.md)。镜像只提供工具链，
不包含流水线源码、配置、星表或观测数据。
