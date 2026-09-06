# Fourier_Quad C++ Pipeline

> English: [README.md](README.md)

## 概述

面向 DECam 数据、基于 MPI 并行的 Fourier_Quad 弱透镜剪切测量 C++17 Pipeline。
统一可执行文件 `Fourier_Quad_Pipe` 提供六个顶层阶段，涵盖星表准备、归档初始化、
剪切星表生成和场畸变（FD）测试。

## Pipeline 版本

| Pipeline | 源码目录 | 用途 |
|---|---|---|
| C++ Lite | [`cpp_Lite`](cpp_Lite/) | 删除未使用替代分支的固定生产路径。 |
| C++ Standard | [`cpp_Standard`](cpp_Standard/) | 保留可选科学路径的完整分支集合。 |

两个版本使用相同的可执行文件名、INI 格式和命令行接口。Lite 从源码中物理删除了
Standard-only 分支，无法通过配置值重新启用。

## 快速开始

### 1. 下载 Release 源码包

打开 [GitHub Releases](https://github.com/Syoong-s/Fourier_Quad_Pipeline_Neo/releases)，
下载准备运行的版本对应源码包：

| Pipeline | Release 源码包 |
|---|---|
| C++ Lite | `cpp-lite-<tag>.zip` |
| C++ Standard | `cpp-standard-<tag>.zip` |

推荐使用固定 Release，以便记录并复现实验实际使用的软件版本。

### 2. 准备运行环境

安装下文[运行环境前置条件](#运行环境前置条件)中的依赖。

### 3. 准备输入数据

准备[输入数据要求](#输入数据要求)中的四类数据。Science images、Gaia catalog 和
External source catalog 为必选；是否需要 DQ masks 取决于所选版本和配置。

### 4. 配置 Pipeline

进入解压后的源码目录，复制完整运行时模板，再编辑路径、数据集、处理阶段和运行时
科学选项。

可用的配置方式以及优先级为：

```text
config文件默认值 < --config INI < 命令行选项CLI
```

未映射到 INI 的固定设置仍位于 `config/*.hpp`，修改后需要重新编译。完整运行时和
编译时参数参考见 [CPP_PIPELINE_PARAMETERS.md](CPP_PIPELINE_PARAMETERS.md)。

*`cpp_Lite` 应先在 `Initialize.hpp` 中设置常用编译默认值；单次运行的覆盖项使用
INI 或 CLI。*

### 5. 编译和运行

```bash
make -j4
./Fourier_Quad_Pipe --help
mpirun -np 4 ./Fourier_Quad_Pipe --config pipeline.ini
```

Docker 与 HPC 运行说明入口：

| 环境 | 指南 |
|---|---|
| Docker 本地运行 | [C++ Docker 指南](cpp_docker/README-CN.md) |
| HPC / Slurm + Apptainer | [C++ runner 指南](cpp_docker/runner/README-CN.md) |

完整编译变量、阶段选择、输入输出、运行方式和失败规则见
[C++ 指南](CPP_GUIDE_CN.md)。

## 运行环境前置条件

| 环境 | 前置条件 | 参考软件栈 |
|---|---|---|
| 常规 Linux | 64 位 Linux；支持 C++17 的 MPI C++ 编译器；OpenMPI 或兼容的 MPI；CFITSIO；FFTW3（双精度和单精度）；Eigen3；BLAS / LAPACK | GCC 12.3.0；OpenMPI 4.1.8；CFITSIO 4.6.4；FFTW 3.3.11；Eigen 3.4.0 |
| HPC / Slurm | 常规 Linux 环境要求；共享文件系统；Slurm；兼容 PMI2 的启动方式；Apptainer / Singularity | 与已发布的 x86_64 容器镜像相同 |

## 输入数据要求

| 输入 | 用途 | 最低要求 |
|---|---|---|
| Science images | 用于源检测、形状测量和弱透镜处理的科学曝光图像。 | 代码支持的 Science FITS/FZ 数据，且符合配置的归档、曝光与 CCD 命名规则。 |
| Gaia catalog | 为天体匹配和测天标定提供天球坐标参考。 | 覆盖 Science images 天区，并符合配置的 Type 1 或 Type 2 分片布局。 |
| External source catalog | 为匹配和下游处理提供天球位置、`zp` 与光度。 | 规范的一度分片，包含 `ra`、`dec`、`zp` 和至少一个观测波段的星等。 |
| DQ masks（*取决于配置*） | 标记坏像素、饱和、探测器缺陷及其他无效像素。 | Lite 以及所有会读取 DQ 数据的 Standard 配置必须提供。 |

### Science images

Science images 是主要科学曝光图像，不是标定星表或输出目录。它们必须采用所选版本
支持的 FITS/FZ 格式，能够被其 FITS 逻辑读取，并符合详细指南中的曝光/CCD 命名和
列表约定。

### Gaia catalog

Gaia catalog 为天体匹配和测天标定提供精确 RA/Dec 参考位置。它必须覆盖实际 Science
images 天区，并直接存放在配置的 `ASTROMETRY_CAT` 目录下。每个被消费的瓦片第一行是
表头，后续每行的前两个数值字段作为 RA 和 Dec；数据行可用逗号或空白分隔，额外字段
会被忽略。

`LensingConfig::AstroCatType` 选择两种文件布局：类型 `1` 是下述旧式大瓦片；类型 `2`
是一度瓦片，可由 `process_astrocat` 从平铺目录中的原始文件生成。

**文件名规范：**
-  `|Dec| < 80°` : `gaia_<p|m><D>_<RR>.cat`
其中`D = floor(|Dec| / 10) + 1`（1-8），`RR = floor(RA / 10)`（00-35，不足两位
补零）。
-  `|Dec| >= 80°` : 不带 RA 后缀 `gaia_<p|m>9.cat`。
- 其中, `p`表示 Dec 非负，`m` 表示 Dec 为负。
- 单个文件必须包含一行表头
> 示例：
> 1. gaia_p1_00.cat 覆盖 `0° <= RA < 10°`、`0° <= Dec < 10°`
> 2. gaia_m3_12.cat 覆盖 `120° <= RA < 130°`、`-30° < Dec <= -20°`
> 3. gaia_p9.cat 覆盖 `0° <= RA < 360°`、`80° <= Dec <= 90°`

*对于类型 1，为了保证位于10°×10°格点边缘的曝光仍能选到星，建议单个星表的dec覆盖上下限在上述基础上
增减2°，ra上下限在0°/30°/60°范围分别增减2°/4°/6°。*

类型 2 文件名为
`<ASTROMETRY_TILE_PREFIX>RA_<RA0>_<RA1>_Dec_<Dec0>_<Dec1>.dat`。
`config/pathconfig.hpp` 中的默认前缀是 `astra_`，前缀本身不含 `RA_`。RA 边界固定
三位，Dec 边界使用 `p`/`m` 加两位绝对值。例如
`astra_RA_123_124_Dec_m05_m04.dat` 覆盖
`123° <= RA < 124°`、`-5° <= Dec < -4°`。

`process_astrocat` 只读取 `--astrocat-input` 的直接普通文件，生成上述类型 2 瓦片，并
删除两个坐标上精确或相差不超过 1 ULP 的重复记录。`--astrocat-output` 只控制发布目录；
如需主流程消费结果，应另行把 `ASTROMETRY_CAT` 指向该目录并设置 `AstroCatType=2`。

### External source catalog

最低 schema 不绑定固定 survey 或波段：

| 字段 | 含义 |
|---|---|
| `ra` | Right Ascension。 |
| `dec` | Declination。 |
| `zp` | 所选 Pipeline 配置实际使用的 catalog `zp` 量。 |
| 任一观测波段星等 | 所选分析使用的任意一个波段 magnitude。 |

用户可以保留额外颜色、redshift、object class、shape 和 flag，但它们不属于最低输入
要求。实际列位置、列名、delimiter、header 处理和投影由 `config/ExtCatConfig.hpp`
或其已记录的 CLI 覆盖项配置。

**文件名规范：**
- 1° × 1° 分片：`<SOURCE_CAT_TILE_PREFIX>RA_<RA0>_<RA1>_Dec_<Dec0>_<Dec1>.dat`
- `config/pathconfig.hpp` 中的默认 `SOURCE_CAT_TILE_PREFIX` 是 `extern_`，前缀本身不含 `RA_`。
- RA 边界固定为三位数字；Dec 边界由 `p` 或 `m` 加两位绝对值组成；每个上边界比下边界大 1 度。
- 单个文件必须包含一行表头
> 示例：
> `extern_RA_123_124_Dec_m05_m04.dat` 覆盖 `123° <= RA < 124°`、`-5° <= Dec < -4°`。

原始下载结果应通过 C++ `process_extcat` 转换为上述流水线规范瓦片；下载器自身的原始
文件名不构成另一套 C++ 命名约定。

### DQ masks（可选）

DQ masks 标记不应参与科学测量的坏像素、饱和、探测器缺陷及其他无效区域。只有当
所选配置关闭 DQ 访问时才可省略。Standard 用户可用 `include_Mask` 选择 mask 模式；
当前 C++ Lite 固定使用逐 CCD DQ masks，因此 Lite 运行仍须提供它们。如果省略 DQ
masks，必须确认所有有效分支和配置路径都不再读取它们。

精确星表文件名、schema、生产端/消费端路径关系和归档约定见
[C++ 指南](CPP_GUIDE_CN.md)。运行时值在 `pipeline.ini` 中配置；
[CPP_PIPELINE_PARAMETERS.md](CPP_PIPELINE_PARAMETERS.md) 标明仍为编译时设置的参数。

## 详细指南

- [C++ Pipeline 指南](CPP_GUIDE_CN.md) /
  [English](CPP_GUIDE.md)
- [C++ 参数参考](CPP_PIPELINE_PARAMETERS.md)
- [C++ Docker 指南](cpp_docker/README-CN.md) /
  [English](cpp_docker/README.md)
- [C++ Slurm/Apptainer runner](cpp_docker/runner/README-CN.md) /
  [English](cpp_docker/runner/README.md)


## 辅助插件

此仓库提供了一份供 Codex/Claude Code 使用的插件，当它们加载插件中技能后，
能够获得下面的能力：

1. 指导使用者配置环境。
2. 提供运行程序的方式。
3. 根据需要设置运行参数。

详情请见 [Codex/Claude Code Plugin](https://github.com/Syoong-s/FQNeoAIManual).


## 许可证

仓库自有代码采用 [MIT License](LICENSE)。容器依赖保留其上游许可，详见
[第三方声明](cpp_docker/THIRD_PARTY_NOTICES.md)。
