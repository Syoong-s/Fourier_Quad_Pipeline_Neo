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
科学选项：

```bash
cp pipeline.example.ini pipeline.ini
```

配置优先级为：

```text
编译默认值 < --config INI < 命令行选项
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

## 许可证

仓库自有代码采用 [MIT License](LICENSE)。容器依赖保留其上游许可，详见
[第三方声明](cpp_docker/THIRD_PARTY_NOTICES.md)。
