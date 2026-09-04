# Fourier_Quad C++ Pipeline

> English: [README.md](README.md)

基于 MPI 并行的 Fourier_Quad 弱透镜 C++17 Pipeline，面向 DECam 数据。统一可执行文件
`Fourier_Quad_Pipe` 可以重分块外部源星表、初始化压缩的 Science/DQ 归档、运行九个
数值阶段、重排结果星表并执行场畸变（FD）剪切测试。

## 选择版本

| 版本 | 适用场景 |
|---|---|
| [`cpp_Standard`](cpp_Standard/) | 需要可选的 flat、mask、astrometry、外部 PSF、混合 PSF 或 PCA PSF 分支。 |
| [`cpp_Lite`](cpp_Lite/) | 使用固定生产路径：Gaia 测天、DQ masks、外部源星表、局域多项式 PSF，且不使用 PCA。 |

两个版本使用相同的可执行文件名和命令行接口。Lite 从源码中删除了不用的分支，并非只
是采用不同默认值的 Standard。

## 快速开始

### 1. 下载 Release 源码包

打开 [GitHub Releases](https://github.com/Syoong-s/Fourier_Quad_Pipeline_Neo/releases)，
只下载准备运行的版本对应的源码包：

| 版本 | Release 源码包 |
|---|---|
| C++ Lite | `cpp-lite-<tag>.zip` |
| C++ Standard | `cpp-standard-<tag>.zip` |

建议使用固定 Release，以便记录并复现分析实际使用的软件版本。

### 2. 准备运行环境

安装下方[运行前置条件](#运行前置条件)列出的依赖。

### 3. 准备输入数据

准备[输入数据要求](#输入数据要求)中的四类输入。Science images、Gaia catalog 和
External source catalog 为必选；是否需要 DQ masks 取决于所选版本和配置。

### 4. 配置 Pipeline

进入解压后的 `cpp_Lite` 或 `cpp_Standard` 目录，复制完整运行时模板，再编辑路径、
数据集、处理阶段和运行时科学选项：

```bash
cp pipeline.example.ini pipeline.ini
```

配置优先级为：

```text
编译默认值 < --config INI < 命令行选项
```

未映射到 INI 的固定数值设置仍位于 `config/*.hpp`，修改后需要重新编译。完整的运行时
与编译时参数参考见 [CPP_PIPELINE_PARAMETERS.md](CPP_PIPELINE_PARAMETERS.md)。

### 5. 编译并运行

```bash
make -j4
./Fourier_Quad_Pipe --help
mpirun -np 4 ./Fourier_Quad_Pipe --config pipeline.ini
```

Docker 与 HPC 运行说明入口：

| 实现 | Docker 本地运行 | HPC / Slurm + Apptainer |
|---|---|---|
| C++ | [C++ Docker 指南](cpp_docker/README-CN.md) | [C++ runner 指南](cpp_docker/runner/README-CN.md) |

完整编译变量、阶段选择、运行方式、输入输出和失败规则见
[CPP_GUIDE_CN.md](CPP_GUIDE_CN.md)。

## 运行前置条件

| 类别 | 前置条件 | 普通 Linux 运行 | HPC / Slurm 运行 | 说明 |
|---|---|---:|---:|---|
| 操作系统 | 64-bit Linux | 必需 | 必需 | 推荐现代 Linux 发行版。 |
| C++ 编译器 | 支持 C++17 的 MPI C++ wrapper | 必需 | 必需 | 固定容器使用 GCC 12.3.0。 |
| MPI | OpenMPI 或兼容的 MPI 实现 | 必需 | 必需 | 固定容器使用 OpenMPI 4.1.8。 |
| CFITSIO | CFITSIO development library | 必需 | 必需 | FITS/FZ 读写；容器版本为 4.6.4。 |
| FFTW3 | 双精度和单精度 FFTW3 库 | 必需 | 必需 | Fourier 变换；容器版本为 3.3.11。 |
| Eigen3 | Eigen3 headers | 必需 | 必需 | C++ 数值运算；容器版本为 3.4.0。 |
| BLAS / LAPACK | BLAS 与 LAPACK | 必需 | 必需 | 容器使用 LAPACK 3.11.0 与 OpenBLAS 0.3.33。 |
| 编译工具 | `make`、shell、标准 GNU 工具 | 必需 | 必需 | 编译和辅助脚本。 |
| 共享文件系统 | 所有 rank/节点可见相同输入输出 | 否 | 必需 | 多节点任务必需。 |
| Slurm | 站点调度器与 `srun` | 否 | Slurm 任务必需 | 批处理分配与启动。 |
| PMI2 兼容启动 | Slurm `pmi2` 支持及兼容 MPI | 否 | 仓库 HPC 容器路径必需 | 用于 Slurm 直接启动 MPI。 |
| Apptainer / Singularity | 无 root 容器运行时 | 否 | 使用仓库 HPC 容器路径时必需 | 计算节点必须可用。 |
| 可写处理目录 | 所有 rank 可写的共享输出/工作目录 | 必需 | 必需 | 源数据与输出不应混放。 |

## 输入数据要求

| 输入 | 必选 | 用途 | 最低要求 |
|---|---:|---|---|
| Science images | 是 | Pipeline 对其执行源检测、形状测量和弱透镜处理的科学曝光图像。 | 代码支持的 Science FITS/FZ 数据，且文件组织符合配置的曝光与 CCD 识别规则。 |
| Gaia catalog | 是 | 为天体匹配和测天标定提供高精度天球坐标参考。 | 覆盖 Science images 天区；每个文件第一行为表头，后续行前两个字段为数值型 `ra`、`dec`。 |
| External source catalog | 是 | 为外部源匹配和下游处理提供天球位置、`zp` 与光度。 | 包含 `ra`、`dec`、`zp` 和至少一个观测波段的星等。 |
| DQ masks | 取决于配置（输入类别可选） | 标记坏像素、饱和、探测器缺陷及其他无效像素。 | 仅当所选版本和配置不读取 DQ masks 时才可省略。 |

### Science images

Science images 是主要科学曝光图像，不是标定星表或输出目录。它们必须采用所选版本
支持的 FITS/FZ 格式，能够被其 FITS 逻辑读取，并符合详细指南中的曝光/CCD 命名和
列表约定。

### Gaia catalog

Gaia catalog 为天体匹配和测天标定提供精确 RA/Dec 参考位置。它必须覆盖实际 Science
images 天区，并直接存放在 `[lensing].astrometry_cat` 目录下（编译默认值为
`ASTROMETRY_CAT`）。Pipeline 消费的每个分片都必须有一行表头，后续每行的前两个
数值字段为 RA 和 Dec；数据行可用逗号或空白分隔，额外字段会被忽略。
`[lensing].astrometry_cat_type` 用于选择以下两种布局：

下述查找公式中的 RA 和 Dec 是从当前 Science CCD 头读取的 WCS 参考坐标。

**Type 1（旧式大分片，默认）：**

- `|Dec| < 80°`：`gaia_<p|m><D>_<RR>.cat`，其中
  `D = floor(|Dec| / 10) + 1`（1-8），`RR = floor(RA / 10)`（00-35，不足两位补零）。
- `|Dec| >= 80°`：不带 RA 后缀的 `gaia_<p|m>9.cat`。
- `p` 表示 Dec 非负，`m` 表示 Dec 为负。
- 单个文件必须包含一行表头。

> 示例：
>
> 1. `gaia_p1_00.cat` 覆盖 `0° <= RA < 10°`、`0° <= Dec < 10°`。
> 2. `gaia_m3_12.cat` 覆盖 `120° <= RA < 130°`、`-30° < Dec <= -20°`。
> 3. `gaia_p9.cat` 覆盖 `0° <= RA < 360°`、`80° <= Dec <= 90°`。

*为了保证位于 10° × 10° 格点边缘的曝光仍能选到星，建议单个星表的 Dec 覆盖上下限
在上述基础上增减 2°；在从绝对 Dec 0°、30°、60° 开始的范围内，RA 覆盖上下限分别
增减 2°、4°、6°。*

**Type 2（1 度分片）：**

- `<ASTROMETRY_TILE_PREFIX>RA_<RA0>_<RA1>_Dec_<Dec0>_<Dec1>.dat`。
  `config/pathconfig.hpp` 中默认前缀为 `astra_`，前缀本身不含 `RA_`。RA 边界固定为
  三位数字，Dec 边界由 `p` 或 `m` 加两位绝对值组成。
- 示例：`astra_RA_123_124_Dec_m05_m04.dat` 覆盖
  `123° <= RA < 124°`、`-5° <= Dec < -4°`。

可选的 `process_astrocat` 阶段会读取原始 Gaia 目录下的直接常规文件，转换为 Type 2
分片；RA、Dec 两列均完全相同或相差不超过一个 ULP 的坐标对会被去重（包括分片
边界），输出包含必需的 `RA    DEC` 表头。其 `[astrocat].output_directory` 或
`--astrocat-output` 只控制生产端输出目录，不会与 `[lensing].astrometry_cat` 校验，
也不会更新后者。后续消费这些结果时，须另行把 `[lensing].astrometry_cat` 指向该
目录，并设置 `[lensing].astrometry_cat_type = 2`。

### External source catalog

最低 schema 不绑定固定 survey 或波段：

| 字段 | 含义 |
|---|---|
| `ra` | Right Ascension。 |
| `dec` | Declination。 |
| `zp` | 所选 Pipeline 配置实际使用的 photometric-redshift 量。 |
| 任一观测波段星等 | 所选分析使用的任意一个波段 magnitude。 |

用户可以保留额外颜色、redshift、object class、shape 和 flag，但它们不属于最低输入
要求。实际列位置、delimiter、header 处理和投影由 INI 的 `[extcat]` 段、对应命令行
选项或 `config/ExtCatConfig.hpp` 中的编译默认值配置。

**文件名规范：**

- 1° × 1° 分片：`<SOURCE_CAT_TILE_PREFIX>RA_<RA0>_<RA1>_Dec_<Dec0>_<Dec1>.dat`。
- `config/pathconfig.hpp` 中默认 `SOURCE_CAT_TILE_PREFIX` 为 `extern_`，前缀本身不含 `RA_`。
- RA 边界固定为三位数字；Dec 边界由 `p` 或 `m` 加两位绝对值组成；每个上边界比
  下边界大 1 度。
- 单个文件必须包含一行表头。

> 示例：`extern_RA_123_124_Dec_m05_m04.dat` 覆盖
> `123° <= RA < 124°`、`-5° <= Dec < -4°`。

### DQ masks（可选）

DQ masks 标记不应参与科学测量的坏像素、饱和、探测器缺陷及其他无效区域。只有当
所选版本和配置关闭 DQ 访问时才可省略。Standard 用户可用
`[lensing].include_mask` 选择 mask 模式；Lite 固定使用逐 CCD DQ masks，因此 Lite
运行仍须提供它们。如果省略 DQ
masks，必须确认所有有效分支和配置路径都不再读取它们。

## Pipeline

六个顶层阶段固定按以下顺序运行：

```text
process_astrocat -> process_extcat -> process_init -> process_main -> process_rearr -> process_fd
```

`process_astrocat` 生成去重的一度 Gaia 瓦片；其输出目录与后续通过
`[lensing].astrometry_cat` 消费的 Gaia 目录相互独立。

`process_main` 包含从预处理到星表合并的九个质数门控阶段。完整阶段和数据约定见
[CPP_GUIDE_CN.md](CPP_GUIDE_CN.md)；运行时和编译时设置分别列于
[CPP_PIPELINE_PARAMETERS.md](CPP_PIPELINE_PARAMETERS.md)。

主要产物包括自动生成的曝光表、逐曝光 `*_all.cat` 星表、重排后的 `subcat_*.cat`
星表和 `fdout/FD_test_comb.dat`。输入归档和星表均在原位置只读访问。

## 容器与 HPC

[`cpp_docker`](cpp_docker/) 提供 x86_64 Docker 工具链镜像。镜像包含编译器和依赖库，
不包含 Pipeline 源码或观测数据；两者均从宿主机挂载。Slurm 环境在确认站点提供 PMI2
启动插件后，使用生产级 [Apptainer runner](cpp_docker/runner/README-CN.md)。

## 文档

| 文档 | 用途 |
|---|---|
| [C++ 指南](CPP_GUIDE_CN.md) / [English](CPP_GUIDE.md) | 编译、配置、运行、输入、输出和失败规则 |
| [参数参考](CPP_PIPELINE_PARAMETERS.md) | 每个配置头文件一张完整 Standard/Lite 对照表，并标明 INI/CLI 覆盖 |
| [容器指南](cpp_docker/README-CN.md) / [English](cpp_docker/README.md) | Docker 镜像与本地容器工作流 |
| [Slurm runner](cpp_docker/runner/README-CN.md) / [English](cpp_docker/runner/README.md) | Apptainer/Singularity 部署 |

## 许可证

仓库自有代码采用 [MIT License](LICENSE)。容器依赖保留其上游许可，详见
[第三方声明](cpp_docker/THIRD_PARTY_NOTICES.md)。
