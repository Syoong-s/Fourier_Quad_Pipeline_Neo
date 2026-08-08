# Fourier_Quad C++ 流水线指南

C++17（`cpp_Standard` / `cpp_Lite`）流水线完整指南：源码结构、流水线阶段、参数、外部星表、编译、运行模式、初始化输出布局、Docker 环境与 HPC runner。项目总览与 Fortran 流水线见
[`README_CN.md`](README_CN.md) 与 [`F77_GUIDE_CN.md`](F77_GUIDE_CN.md)。完整参数参考见
[`CPP_PIPELINE_PARAMETERS.md`](CPP_PIPELINE_PARAMETERS.md)。

`cpp_Standard` 为完整构建（含 PCA `PSFRecons`）；`cpp_Lite` 为冻结分支的精简构建（移除 `PSFRecons`）。精简版变更日志见 `cpp_Lite/REFACTOR_NOTES.md`。

> **English**: see [CPP_GUIDE.md](CPP_GUIDE.md)

---

## 源码结构

### 源码目录

#### `cpp_Standard/` - 完整 C++17 流水线

| 文件 | 说明 |
|---|---|
| `main.cpp` | MPI 入口、工作流参数解析与四个函数的执行顺序。 |
| `include/ProcessConfig.hpp` | 共享工作流默认值与函数开关。 |
| `src/process_init/`、`include/process_init/` | 归档初始化器源码与头文件。 |
| `src/process_main/process_main.cpp`、`include/process_main/process_main.hpp` | 曝光列表读取与阶段 1–9 调度。 |
| `src/process_rearr/`、`include/process_rearr/` | 自包含的 `_all.cat` 天区切分、MPI 重分配、排序子星表与汇总输出。 |
| `include/process_rearr/ProcessRearrConfig.hpp` | 重排专属参数，以及 `外部列数 + 1 + ichi2` 派生行宽。 |
| `include/process_main/LensingConfig.hpp` | 配置常量（等价于 `para.inc` + `cust_para.inc` + `sig_para.inc`）。 |
| `src/process_main/PreProcess.cpp`、`include/process_main/PreProcess.hpp` | **阶段 1**：预处理。 |
| `src/process_main/Astrometry.cpp`、`include/process_main/Astrometry.hpp` | **阶段 2**：天体测量校准。 |
| `src/process_main/SourceExtractor.cpp`、`include/process_main/SourceExtractor.hpp` | **阶段 3**：源检测与提取。 |
| `src/process_main/FourierTransformSt1.cpp`、`include/process_main/FourierTransformSt1.hpp` | **阶段 4**：第一阶段傅里叶变换。 |
| `src/process_main/PSFModel.cpp`、`include/process_main/PSFModel.hpp` | **阶段 5**：PSF 建模。 |
| `src/process_main/PSFRecons.cpp`、`include/process_main/PSFRecons.hpp` | PSF PCA 重建（仅 `PSF_Ms=1` 时启用）。 |
| `src/process_main/FourierTransformSt2.cpp`、`include/process_main/FourierTransformSt2.hpp` | **阶段 6**：第二阶段傅里叶变换。 |
| `src/process_main/ShearMeasurement.cpp`、`include/process_main/ShearMeasurement.hpp` | **阶段 7**：Fourier\_Quad 剪切估计。 |
| `src/process_main/ExposureInfo.cpp`、`include/process_main/ExposureInfo.hpp` | **阶段 8**：单次曝光统计。 |
| `src/process_main/CatalogCombiner.cpp`、`include/process_main/CatalogCombiner.hpp` | **阶段 9**：星表合并与标定。 |
| `src/process_main/` 与 `include/process_main/` 中的支撑模块 | FITS 读写、线性代数、图像处理、MPI 调度与通用数值工具。 |
| `src/process_main/UniversalUtils.cpp`、`include/process_main/UniversalUtils.hpp` | 通用工具函数。 |
| `src/process_main/ImageProcessing.cpp`、`include/process_main/ImageProcessing.hpp` | 图像处理工具。 |
| `src/process_main/NumericalRecipes.cpp`、`include/process_main/NumericalRecipes.hpp` | Numerical Recipes 移植（随机数、排序、插值）。 |
| `src/process_main/MPIScheduler.cpp`、`include/process_main/MPIScheduler.hpp` | MPI 初始化与任务分发。 |
| `src/process_main/ExStar.cpp`、`include/process_main/ExStar.hpp` | 恒星提取与分类。 |
| `Makefile` | 构建文件。使用 `mpicxx`，C++17，链接 CFITSIO、FFTW、LAPACK。 |

#### `cpp_Lite/` - 精简 C++17 流水线

采用与 `cpp_Standard/` 相同的 `process_extcat` / `process_init` /
`process_main` / `process_rearr` 集成目录结构和运行参数接口，但科学模块仍保留
冻结后的 Lite 分支，且不含
`PSFRecons.cpp/.hpp`。详见 `cpp_Lite/REFACTOR_NOTES.md`。



## 源码目录布局

- `include/ProcessConfig.hpp`：工作流默认值与默认阶段开关。
- `include/process_extcat/`、`src/process_extcat/`：外部星表 schema、解析、MPI 字节范围切分与确定性分片发布。
- `include/process_init/`、`src/process_init/`：初始化封装及保留的 `Initializer` 与 `FitsExtractor` 模块。
- `include/process_main/`、`src/process_main/`：`LensingConfig`、全部数值模块、曝光表加载与完整的 Stage 1–9 编排。
- `include/process_rearr/`、`src/process_rearr/`：自包含的 `_all.cat` schema 校验、全天加权 k-d 划分、MPI 重分布、排序子星表发布与汇总输出。
- 每个 `cpp_Standard` / `cpp_Lite` 根目录仅包含可执行入口、构建文件、文档与各阶段实现子树。

## 流水线阶段

流水线共 9 个阶段。阶段执行由 `PROCESS_stage` 参数（定义于 `para.inc` /
`LensingConfig.hpp`）控制，该参数为素因数之积。当 `PROCESS_stage` 能被某阶段的
素数整除时，该阶段执行。默认值 `2 * 3 * 5 * 7 * 11 * 13 * 17 * 19 * 23` 启用全部
阶段。

| 阶段 | 素数 | 函数 | 说明 |
|---|---|---|---|
| 1 | 2 | `pre_process` / `PreProcess` | 读取 FITS 图像，平场与掩膜校正，背景噪声估计（F6 mode-bar 估计器）。 |
| 2 | 3 | `proc_astrometry` / `Astrometry` | 基于 Gaia 参考星表的天体测量校准，WCS 拟合。 |
| 3 | 5 | `proc_source` / `SourceExtractor` | 源检测、去混叠、stamp 提取。 |
| 4 | 7 | `proc_FFT_st1` / `FourierTransformSt1` | 星系 stamp 第一阶段傅里叶变换。 |
| 5 | 11 | `proc_PSF` / `PSFModel` | 基于恒星 stamp 的 PSF 建模（局域多项式拟合）。可选 PCA 重建（`PSF_Ms=1`）。 |
| 6 | 13 | `proc_FFT_st2` / `FourierTransformSt2` | 第二阶段傅里叶变换。 |
| 7 | 17 | `proc_shear` / `ShearMeasurement` | 基于四阶傅里叶矩的 Fourier\_Quad 剪切估计。 |
| 8 | 19 | `proc_info` / `ExposureInfo` | 收集单次曝光统计（PSF FWHM、恒星数等）。 |
| 9 | 23 | `proc_combine_shear_catalog` / `CatalogCombiner` | 跨曝光合并剪切星表并应用标定校正。 |

如需跳过某阶段，将 `PROCESS_stage` 除以该阶段的素因数即可。例如设置
`PROCESS_stage = 2 * 3 * 5 * 7 * 11 * 13 * 17 * 19`（省去 23）可跳过星表合并。

---


## 参数配置

### C++（`LensingConfig.hpp`）

`cpp_Standard` 将三个 Fortran include 文件中的全部参数整合到
`include/process_main/LensingConfig.hpp`。`cpp_Lite` 使用相同的相对路径，但只保留
冻结后的参数子集，详见 `cpp_Lite/REFACTOR_NOTES.md`。星表列索引已调整为 C++ 的
0 基。全部运行参数、外部输入、`ProcessConfig.hpp` 默认值、
`LensingConfig.hpp` 参数、可选值及 Standard/Lite 差异见
[`C++ Pipeline External Inputs and Parameter Reference`](CPP_PIPELINE_PARAMETERS.md)。

---



## 外部源星表

当 `ext_cat = 1` 时，C++ 和 Fortran 流水线会从 `SOURCE_CAT` 指定的
目录读取 1° 网格的 DES Y6 GOLD 星表。
[`gen_src_cat`](gen_src_cat/README.md) 工具可下载或重新分区这些数据，生成
符合流水线文件名约定的分块文件。Python 下载器固定写出 DES Y6 GOLD 的
18 列格式；C++ MPI 重分区器默认保留原始星表的全部列，也可在开启显式列选择后
按任意给定顺序抽取字段。该重分区器也已作为 `process_extcat` 集成进
`cpp_Standard` 与 `cpp_Lite`，并作为 `process_init`、`process_main` 之前的
可选第一阶段。

流水线生成的曝光级 `_all.cat` 还可由自包含的第四个函数 `process_rearr` 按天区
重新切分。直通模式的外部列宽来自 `EXTCAT_TOTAL_COLUMNS`（默认 18）；专属配置头
按 `18 + 1 + ichi2(25) = 44` 派生默认总列数。关闭显式投影时直接使用配置的一基
RA/Dec 原始列号，开启后会按投影列表自动换算位置。输出默认位于每个数据集的
`rearranged_catalog/` 目录。

对于 C++ 外部星表，请在所选流水线的 `ProcessConfig.hpp` 中设置一基原始列号
`EXTCAT_RA_COLUMN_ONE_BASED`、`EXTCAT_DEC_COLUMN_ONE_BASED` 和
`EXTCAT_ZP_COLUMN_ONE_BASED`。默认值分别为 `5`、`6`、`17`，对应 DES Y6 GOLD
中的 `ra`、`dec`、`dnf_z`。`process_main` 只将这三列转换为数值，其他字段可为
任意字符串，也不再要求固定 18 列格式。开启显式投影时必须选入这三列，读取器会
根据投影列表顺序自动换算它们在输出星表中的位置。运行时也可分别使用
`--extcat-ra-column`、`--extcat-dec-column`、`--extcat-zp-column` 覆盖。
`process_rearr` 使用相同的 RA/Dec 规则，但要求完整 `_all.cat` 行的所有字段均为
有限数值。

```bash
cd gen_src_cat
python3 -m venv .venv
source .venv/bin/activate
python -m pip install numpy pyvo
python query_y6gold_sync_mp_v2.py
```

运行前请在脚本中检查天区范围、行数上限、输出目录和并发数。C++ 程序应在
`LensingConfig.hpp` 中设置主路径 `SOURCE_CAT`，也可通过
`--extcat-output` 覆盖；`EXTCAT_OUTPUT_DIRECTORY` 会跟随该路径，使生成阶段与
数值流水线读取同一目录。已有的任意数量原始星表也可通过
`--run-extcat true --extcat-input PATH` 在流水线内并行分块。Fortran 仍在
`para.inc` 中设置 `SOURCE_CAT`。Python 列格式、
C++ 投影模式、并行行为及下载脚本限流说明详见
[`gen_src_cat/README.md`](gen_src_cat/README.md)。

---



## 编译器与依赖库

本地验证使用 GCC/G++ 15.2.0、Open MPI 5.0.10、CFITSIO 4.6.3、FFTW 3.3.10、
Eigen 3.4.0、OpenBLAS/LAPACK 0.3.33（WSL2 下的 Linux）。

可移植 HPC 目标为 GCC/G++ 12.3.0、Open MPI 4.1.8、CFITSIO 4.6.4、FFTW 3.3.11、
Eigen 3.4.0、LAPACK 3.11.0。集群可使用等价的站点 module，只要存在一个可编译并启动程序的 MPI C++ wrapper。

## 编译

当所有头文件与库都在编译器常规搜索路径中时：

```bash
make -j4
```

若使用一个合并的科学栈前缀（头文件、`include/eigen3` 下的 Eigen、`lib` 下的库）：

```bash
make STACK_PREFIX="${STACK_PREFIX}" -j4
```

当 Eigen 与链接库使用不同前缀时，显式传入可移植的本地覆盖：

```bash
make CXX="${MPI_PREFIX}/bin/mpicxx"   STACK_PREFIX="${STACK_PREFIX}"   EIGEN_INCLUDE="${EIGEN_INCLUDE}" -j4
```

无需 Windows 原生编译器或 wrapper。集群构建在加载站点的编译器、MPI 与科学库 module 后使用同一 Makefile。

运行外部星表列读取器测试：

```bash
make test-extcat-reader
```

运行自包含的重排单元与 MPI 集成测试：

```bash
make test-rearr
```

## 默认值与选项语法

编辑 `include/ProcessConfig.hpp` 设置常规数据集与默认执行模式。`RUN_PROCESS_EXTCAT`、`RUN_PROCESS_INIT`、`RUN_PROCESS_REARR` 默认为 `false`；`RUN_PROCESS_MAIN` 默认为 `true`。每个命令行选项都可省略，并覆盖其配置默认值。`--name value` 与 `--name=value` 均可接受，顺序任意。

`EXTCAT_*` 值配置原始星表发现、输出目录、解析策略、MPI 任务大小、可选有序列选择与原始 RA/Dec/ZP 列。禁用显式选择时，每个原始输入字段原位保留。`process_main` 只读取这三个配置字段；未选中的星表列无需为数值。`EXTCAT_TOTAL_COLUMNS` 给 `process_rearr` 提供透传外部宽度；显式投影则使用投影列表长度（即发射的外部 schema）。主星表路径在 `include/process_main/LensingConfig.hpp` 中用 `SOURCE_CAT` 设置。`EXTCAT_OUTPUT_DIRECTORY` 是该值的只读引用，故 `process_extcat` 写入处即 `process_main` 读取处。`--extcat-output` 对单次调用覆盖两者。

`DATASETS` 存储成对的 target/prefix，`CONTAINS` 存储按 OR 语义匹配的归档 basename token。例如：

```cpp
inline const std::vector<DatasetSpec> DATASETS = {
    {"g2013", "c4d_13"},
    {"g2014", "c4d_14"},
    {"g2019", "c4d_19"},
};
inline const std::vector<std::string> CONTAINS = {"v1", "v2"};
```

| 选项 | 用途 |
|:---|:---|
| `--run-extcat BOOL` | 启用/禁用外部星表重分区。 |
| `--run-init BOOL` | 启用/禁用运行时归档初始化。 |
| `--run-main BOOL` | 启用/禁用运行时数值流水线。 |
| `--run-rearr BOOL` | 启用/禁用 `_all.cat` 空间重排；与 `process_main` 同时启用时在其后运行。 |
| `--extcat-input PATH` | 包含任意数量原始文本星表的目录。 |
| `--extcat-output PATH` | 目标分片目录，即有效的 C++ `SOURCE_CAT`。 |
| `--extcat-contains TEXT` | 可重复、区分大小写的 basename token；任一匹配即选中。 |
| `--science-root PATH` | 只读多 HDU Science FITS/FZ 归档库。 |
| `--dq-root PATH` | 与 Science 配对的只读 DQ FITS/FZ 归档库。 |
| `--output-root PATH` | 目标目录与生成 list 的父目录。 |
| `--dataset TARGET:PREFIX` | 可重复的 target/prefix 对；首个 `--dataset` 清除编译期列表。 |
| `--contains TOKEN` | 可重复的 basename token（OR 匹配）。 |
| `--existing fail\|resume\|overwrite` | 既有输出策略。 |
| `--f77-max-path N` | 仅 initializer 使用的生成路径长度上限；`0` 禁用检查。 |
| `--expo-list PATH` | 单条曝光表（main/rearr-only 模式）。 |

`--f77-max-path` 只属于 `process_init`：默认值 150 用于保持与旧 Fortran 工作流的路径兼容性，
不会截断或拒绝 `process_main` 的路径。main 中的路径是 `std::string`，其实际边界来自所用
文件系统和 I/O 库（FITS 产物还受 CFITSIO 文件名容量约束）。

## 运行模式

仅外部星表执行接受任意数量的匹配输入文件，且不要求配置数据集：

```bash
mpirun -np 4 ./Fourier_Quad_Pipe   --run-extcat true --run-init false --run-main false   --extcat-input /data/raw_catalogs   --extcat-output /data/catalogs/des_y6_chunks   --extcat-contains .csv --extcat-contains y6_gold
```

`process_extcat` 在数据集循环前集体运行一次。其默认输出保留完整原始 schema；例如 `--extcat-columns 5,3,4,1` 将原始 5、3、4、1 列写为输出 1–4 列。若失败，后续阶段不启动。传递给 `process_main` 的输出必须包含 `EXTCAT_RA_COLUMN_ONE_BASED`、`EXTCAT_DEC_COLUMN_ONE_BASED`、`EXTCAT_ZP_COLUMN_ONE_BASED` 配置的原始列；仅重排运行需要 RA 与 Dec，无需 ZP。

仅 main 本地执行：

```bash
mpirun -np 4 ./Fourier_Quad_Pipe   --run-init false --run-main true   --expo-list /data/work/expo_g2019.list
```

仅重排本地执行使用同一曝光表引用的既有 per-exposure `_all.cat`：

```bash
mpirun -np 4 ./Fourier_Quad_Pipe   --run-init false --run-main false --run-rearr true   --expo-list /data/work/expo_g2019.list
```

仅初始化本地执行：

```bash
mpirun -np 4 ./Fourier_Quad_Pipe   --run-init true --run-main false   --science-root /data/archive/science   --dq-root /data/archive/dq   --output-root /data/work --dataset g2019:c4d_19
```

批量初始化将每个 target 与其 prefix 配对。多个 contains token 在任一 token 出现于 basename 时选中该归档：

```bash
mpirun -np 4 ./Fourier_Quad_Pipe   --run-init true --run-main false   --science-root /data/archive/science   --dq-root /data/archive/dq   --output-root /data/work   --dataset g2013:c4d_13 --dataset g2014:c4d_14   --dataset g2019:c4d_19   --contains v1 --contains v2
```

链式本地执行使用相同初始化选项并启用下游阶段开关。初始化成功后，`process_main` 与 `process_rearr` 接收 `process_init` 返回的归一化绝对路径 `output_root/expo_<target>.list`。该生成路径覆盖 `--expo-list`、遗留位置参数及所有配置的曝光表默认值。

```bash
mpirun -np 4 ./Fourier_Quad_Pipe   --run-init true --run-main true --run-rearr true   --science-root /data/archive/science   --dq-root /data/archive/dq   --output-root /data/work --dataset g2019:c4d_19   --existing resume
```

启用全部四个开关可先构建外部星表，再初始化、处理并重排每个配置数据集。有效的 `--extcat-output` 路径也用于数值源提取器。

数据集在同一 MPI 通信器上顺序执行，并在首个失败处停止。在 main/rearr-only 批量模式下省略 `--expo-list`；驱动为每个数据集派生一个 `output_root/expo_<target>.list` 路径。单条外部曝光表仅可用于单个下游-only 数据集。在链式批量模式下，每个初始化生成的绝对 list 覆盖该数据集的外部曝光表输入。

在 Slurm 集群上，用站点启动器使用相同的可执行参数，例如 `srun -n 40 ./Fourier_Quad_Pipe ...`。初始化失败状态是集体的；初始化失败后永不进入数值阶段。

## 初始化输出契约

对每个 `--output-root OUTPUT --dataset TARGET:PREFIX`，初始化严格按以下顺序执行：幂等创建固定类型目录；抽取 Science/DQ chip；为每个成功的 Science 曝光写 `stamps/<EXPOSURE>.list`；由 rank 0 发布两个顶层 list；从已发布的 expo list 创建 chip 产品曝光子目录；最后发布 manifest。源 `.fits.fz` 归档原位读取，从不复制或删除。两个变体均在 `include/OutputLayout.hpp` 中集中声明完整固定目录契约及其 chip 产品子集。

### `OUTPUT/TARGET` 下的输出布局

- `science/<EXPOSURE>/<EXPOSURE>_<N>.fits` - Science chip 图像，按曝光分片，每曝光一个子目录；`<N>` 为二维 HDU 出现序号（1, 2, 3, ...）。
- `dqmask/<EXPOSURE>/<EXPOSURE>_<CCDNUM>.fits` - DQ mask chip 图像，按曝光分片；`<CCDNUM>` 为 FITS `CCDNUM` 头值。
- `stamps/` - 抽取期间写入的 per-exposure chip 表（`<EXPOSURE>.list`），以及 `process_main` 全部中间产物（按类型分子目录，见下）。
- `astrometry/dat_Astro/`、`astrometry/Head/`、`astrometry/dat_Chk/` - 天体测量解（`<P>_astro.dat`）、WCS `.head` 文件与校验数据。
- `result/` - 最终 per-exposure 产物，含 `process_rearr` 消费的 `<EXPOSURE>_all.cat`。

### process_main 中间产物（位于 `stamps/` 下）

按类型分子目录取代了旧的扁平 `stamps/`、`rescale/`、`starxy/`、`fits_psfresi/`、`dat_pcs/`、`dat_starcomp/` 目录：

`Norm/`、`cat_Orig/`、`dat_StarInfo/`、`dat_StarCanInfo/`、`dat_SrcInfo/`、
`dat_PsfFit/`、`dat_Shear/`、`dat_ExpoInfo/`、`dat_StarComp/`、`dat_StarCompV2/`、
`dat_Rescale/`、`dat_StarXY/`、`dat_Pcs/`、`fits_StarCan/`、`fits_StarCanN/`、
`fits_StarCanP/`、`fits_StarP/`、`fits_Src/`、`fits_Noise/`、`fits_SrcP/`、
`fits_PsfLocal/`、`fits_PsfSrc/`、`fits_PsfResi/`。

所有 chip 级产物均在类型目录下再按曝光分一层：
`<类型>/<EXPOSURE>/<CHIP><后缀>`。例如归一化 chip 写入
`stamps/Norm/<EXPOSURE>/<CHIP>_norm.fits`，对应天体测量解写入
`astrometry/dat_Astro/<EXPOSURE>/<CHIP>_astro.dat`。`.head`、
`_star_info_expo.dat`、`_star_power_expo.fits`、`_PSF_source.fits`、
`_expo_info.dat`、`_all.cat` 等曝光级产物仍直接位于原有类型目录中。rank 0 仅在发布并重新读取 `expo_TARGET.list` 后，才幂等创建这些 chip 产品的 `<EXPOSURE>/` 目录。

### process_main 路径与输出失败契约

已按上述布局审计 Stage 1--9 的生产者/消费者链。chip 产物在写入与读取两侧均通过
`OutputLayout::chipPath` 构造；曝光级产物仍直接位于其声明的类型目录。已核对的链路包括：
天体测量/归一化、WCS/check、源与星候选提取、星功率、PSF 产物、源功率、shear、曝光信息和
最终星表合并。DQ 输入读取路径与 initializer 发布的
`dqmask/<EXPOSURE>/<EXPOSURE>_<CCDNUM>.fits` 契约一致；未发现目录层级或后缀不匹配。

`process_main` 的全部文本输出统一使用 checked `MainIO::OutputFile`，FITS 输出统一由 checked
`FitsIO` 创建、写入和关闭。create、write、flush 或 close 任一失败时，程序输出一条
`Output creation failed` 信息，包含 MPI rank、操作、输出路径及 OS/CFITSIO 原因，随后调用
`MPI_Abort` 终止整个作业，避免 master 或其他 worker 留在动态调度等待中。

### 曝光表生成

抽取期间每个 rank 直接从抽取结果写出 `stamps/<EXPOSURE>.list`（其产出的 chip 图像路径）- 不进行事后磁盘反查。所有 rank 完成后，rank 0 扫描 `stamps/`、排序 per-exposure 表，并先原子发布：

- `OUTPUT/expo_TARGET.list` - 顶层表；每行为 `"<OUTPUT/TARGET/stamps/<EXPOSURE>.list>"  <chip 数>`。
- `OUTPUT/fits_TARGET.list` - 全部 Science chip 图像路径的扁平表。

随后 rank 0 读取已发布的 expo list 并创建各 chip 产品曝光子目录；该步骤成功后才原子发布 `OUTPUT/init_TARGET_manifest.json`。清单 schema 版本 2 记录 `exposure_directories_created` 完成标志，并在 `filename_tokens` 数组中记录所有启用的 basename 过滤器。Science chip 按二维 HDU 出现编号；DQ chip 按 `CCDNUM` 编号。下游阶段通过 `getDir(image, 3)`（向上三层：`science/<EXPOSURE>/<file>` -> `science` -> `OUTPUT/TARGET`）从 Science chip 路径推导 dataset root，per-chip DQ mask 从 `dqmask/<EXPOSURE>/<EXPOSURE>_<CCDNUM>.fits` 读取。

## Docker 环境


Docker 环境提供可复现的构建工具链，无需手动安装编译器和科学库。


### Docker 目录

#### `cpp_docker/`

构建可移植 HPC 运行时镜像：

| 组件 | 版本 |
|---|---|
| 基础镜像 | Rocky Linux 8.10 |
| G++ | 12.3.0 (conda-forge) |
| OpenMPI | 4.1.8 (Slurm PMI2 直启动) |
| CFITSIO | 4.6.4 |
| FFTW | 3.3.11 |
| Eigen | 3.4.0 |
| LAPACK / OpenBLAS | 3.11.0 / 0.3.33 |

关键文件：`Dockerfile`、`compose.yaml`、`pixi.toml`、`.env.example`、
`scripts/verify-image.sh`、`scripts/check-public-repo.sh`、`runner/`（HPC 部署
脚本）、`SOURCES.md`、`THIRD_PARTY_NOTICES.md`。

---


### 快速开始（C++）

```bash
cd cpp_docker
cp .env.example .env
# 编辑 .env：设置 CPP_SOURCE_HOST、星表路径和 PROCESS_DATA_HOST

docker compose build
docker compose run --rm FourierQuad-CPP
```

进入容器后：

```bash
cd /workspace/src_pipe
make -j4
mpirun -np 4 ./Fourier_Quad_Pipe /data/DataProcess/expo_list.list
```


### 验证镜像

每个 Docker 目录附带验证脚本：

```bash
bash f77_docker/scripts/verify-image.sh f77pipeline-dev:gnu4.8.5
bash cpp_docker/scripts/verify-image.sh cpppipeline-dev:gxx12.3-openmpi4.1.8-pmi2
```

Docker 环境详细文档请参见：
- [`f77_docker/README.md`](f77_docker/README.md) / [`f77_docker/README-CN.md`](f77_docker/README-CN.md)
- [`cpp_docker/README.md`](cpp_docker/README.md) / [`cpp_docker/README-CN.md`](cpp_docker/README-CN.md)

---



## HPC 部署

两套 Docker 环境均包含 `runner/` 目录，提供 Slurm/Apptainer 部署脚本。典型流程：

1. 拉取 GHCR 镜像并转换为 SIF：
   ```bash
   bash f77_docker/runner/pull-sif.sh    # 或 cpp_docker/runner/pull-sif.sh
   ```

2. 配置环境：
   ```bash
   cp f77_docker/runner/f77pipeline.env.example f77_docker/runner/f77pipeline.env
   # 编辑集群路径
   ```

3. 审计 MPI 兼容性（只读）：
   ```bash
   bash f77_docker/runner/inspect-cluster-mpi.sh
   ```

4. 提交流水线：
   ```bash
   sbatch f77_docker/runner/f77pipeline.slurm
   ```

HPC runner 详细文档请参见：
- [`f77_docker/runner/README.md`](f77_docker/runner/README.md) / [`f77_docker/runner/README-CN.md`](f77_docker/runner/README-CN.md)
- [`cpp_docker/runner/README.md`](cpp_docker/runner/README.md) / [`cpp_docker/runner/README-CN.md`](cpp_docker/runner/README-CN.md)

---
