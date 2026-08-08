# cpp_Standard 单 SIF Docker 与 runner 开发总结

## 1. 目标与边界

本实现沿用 `f77_docker` 的制品边界：

- 镜像只包含编译器、MPI、科学库和运行所需系统文件；
- `cpp_Standard` 源码、压缩 Science/DQ 归档、星表、平场、FITS 数据和输出位于镜像外；
- 本地使用 Docker bind，集群使用一个预编译 Apptainer SIF；
- 编译发生在挂载的源码副本中，产物保留在共享文件系统；
- 真实 pipeline 前必须依次通过镜像、完整编译、单节点 MPI 和多节点 MPI
  验证。

实现不再区分“本地镜像”和“pilogin 镜像”。pilogin 只是验证站点，不是
镜像 profile。支持契约是 x86_64 Linux、Slurm PMI2、Apptainer/Singularity
和同路径共享文件系统。

## 2. 唯一工具链

最终制品固定为：

- G++ 12.3.0；
- OpenMPI 4.1.8；
- Slurm 25.11.2 源码构建的最小 `libpmi2` client；
- CFITSIO 4.6.4；
- FFTW 3.3.11；
- Eigen 3.4.0；
- LAPACK/BLAS 3.11.0；
- OpenBLAS 0.3.33。

OpenMPI 5 不提供项目需要的 PMI2 路径，而许多现有 Slurm 集群仍只暴露
PMI2。这里选用稳定的 OpenMPI 4.1.8，从官方归档构建，并通过
`--with-pmi=/opt/pmi2` 链接镜像内 PMI2 client。OpenMPI configure 已实际
检测到 `pmi2.h`、`libpmi2` 和 `PMI2_Init`；最终镜像的
`ompi_info --all` 包含 `MCA ess: pmi`。OpenMPI 同时使用
`--with-slurm` 启用 `srun` direct-launch 所需的 Slurm 环境识别组件，
但不链接宿主 `libslurm`。

Slurm 源码只用于构建 `contribs/pmi2`。镜像不安装 Slurm controller、
daemon、命令集、配置或 MUNGE。PMI2 client 是应用与宿主 Slurm step 的
协议边界，不要求宿主 OpenMPI 与容器 OpenMPI 同版本。

## 3. Docker 和 Pixi 结构

Dockerfile 只有一个最终 target：`runtime`。Pixi lock 只有一个环境：
`default`，且不包含 Conda OpenMPI 或 PMIx，避免重复 MPI 污染。

构建阶段依次为：

1. 从 digest-pinned Pixi 镜像取得可执行文件；
2. 在 digest-pinned Rocky Linux 8.10 上安装锁定的编译器和科学栈；
3. 编译最小 Slurm PMI2 client；
4. 使用同一 G++/binutils 工具链编译 OpenMPI 4.1.8；
5. 只复制安装结果到非 root runtime。

基础镜像、上游 URL 和 SHA256 记录在 `SOURCES.md`。最终 runtime 不包含
下载归档、解包源码或 pipeline 源码。

## 4. cpp_Standard 编译接口

Makefile 使用标准 `CPPFLAGS`、`CXXFLAGS`、`LDFLAGS` 和 `LDLIBS`，
编译器为 `mpicxx`，语言标准为 C++17。`STACK_PREFIX` 是可选外部前缀；
容器默认依靠自身环境变量，无需写死站点路径。

Makefile 生成 `.d` 自动依赖文件，`clean` 只删除当前源码目录的 object、
dependency 和 `Fourier_Quad_Pipe`。`LensingConfig.hpp` 中的星表和平场
路径属于科学配置，本实现没有修改；runner bind 目标必须与它一致。

## 5. 通用 Slurm runner

集群启动链为：

`srun --mpi=pmi2` → `run-apptainer.sh` → `apptainer exec --cleanenv` →
SIF 内 OpenMPI 应用。

runner 不调用宿主 `mpirun`，也不加载宿主 MPI。`--cleanenv` 清除
`OMPI_*`、`PMIX_*`、`LD_LIBRARY_PATH` 等潜在污染，然后显式转发 Slurm
生成的全部 `SLURM_*`、`PMI_*`、`PMI2_*`，以及配置 allowlist。

作业职责被明确拆分：

- `build-sif.slurm` / `pull-sif.sh`：临时输出、SHA256 sidecar、原子改名；
- `run-apptainer.sh --check`：SIF hash、版本、PMI 组件和 bind；
- `compile-pipeline.slurm`：单 task 编译完整 pipeline；
- `mpi-smoke-test.slurm`：同一模板执行单节点或多节点 PMI2 smoke；
- `cpppipeline.slurm`：只启动已编译的真实 pipeline。

站点只需修改 env 中的绝对共享路径、资源参数、可选 runtime module 和网络
调优；这些不是新的镜像或 runner 分支。

## 6. 可移植性说明

基线 OpenMPI 保留节点内共享内存和节点间 TCP，不绑定 UCX、OFI、
InfiniBand verbs 或 vendor provider。因此一个 SIF 可以避开站点驱动 ABI
差异，优先保证可运行性。高性能网络加速必须由站点单独验证，不能仅凭
TCP smoke 宣称。

一个 x86_64 Slurm/PMI2 SIF不能自然覆盖 ARM、PBS/LSF 或只提供 PMIx 的
站点。这些属于不同运行契约，而不是通过增加 `local`/`pilogin` 镜像解决。

## 7. 本地验证证据

- 最终 Docker manifest-list ID：
  `sha256:3640b5add63167aa6a5919939f5139c6ca56b6288050e1a66935935009501b47`；
- Docker 报告平台为 linux/amd64，镜像大小 376,237,960 bytes；
- Docker archive 为 376,257,536 bytes，SHA256
  `07c0277647a04d92533465a48c34ed9872d8c9a388cd0398b4492c8aef3c0304`；
- 精确 G++ 12.3.0、OpenMPI 4.1.8、`ess:pmi`、`libpmi2`、CFITSIO、FFTW、
  Eigen、LAPACK/BLAS 检查通过，且存在 `MCA schizo: slurm`；
- 两个 Docker-local MPI rank 和科学库综合 smoke 通过；
- 完整 18 个 translation unit 在隔离源码副本中重新编译和链接成功；
- `ldd` 显示 `libmpi.so.40` 来自 `/opt/openmpi-4.1.8`，
  `libpmi2.so.0` 来自 `/opt/pmi2`，科学库来自 Pixi stack；
- 只有原有源码的未使用参数/变量和 Eigen 模板 warning，没有编译或链接
  error；
- Bash 语法与 `scripts/check-public-repo.sh` 通过。

Docker 内部自测的 `mpirun` 使用 OpenMPI 内置 `plm=isolated`，仅用于验证
容器自包含的双 rank；它不属于集群执行路径。

## 8. pilogin 验证位置

pilogin 验证工作区为：

`~/ysx/Agents/cpppipeline-portable`

目录和子目录在写入前完成 `readlink`/owner/mode 检查，均解析到授权
`Agents` 子树。所有 Slurm stdout/stderr 写入：

`~/ysx/Agents/outputs`

最终验证制品：

- SIF 大小 336,166,912 bytes；
- SIF SHA256：
  `0db2ada5eab77d243f49a998043e55a5ba620ac0ad93690817ffb4d7dae1591c`；
- SIF sidecar 通过 `sha256sum --check`。

最终 Slurm 记录：

| 门槛 | Job ID | 状态 | 节点 | 用时 | 结果 |
|---|---:|---|---|---:|---|
| 修正版 SIF 构建 | `60090091` | `COMPLETED 0:0` | `cas008` | 7:38 | SIF 与 sidecar 原子落盘 |
| 完整 pipeline 编译 | `60090603` | `COMPLETED 0:0` | `cas008` | 0:25 | 18 个 translation unit；可执行文件 958,888 bytes |
| 单节点 PMI2 smoke | `60090644` | `COMPLETED 0:0` | `cas008` | 0:10 | 2 ranks；科学栈通过 |
| 双节点 PMI2 smoke | `60090690` | `COMPLETED 0:0` | `cas235`,`cas305` | 0:13 | 4 ranks，按 2+2 分布；科学栈通过 |

pilogin 要求多节点作业使用独占节点或每节点申请 40 tasks，因此最终双节点
作业在提交时增加了 `--exclusive`。这是站点资源参数，不是镜像或 runner
分支；实际 MPI 进程仍为每节点 2 ranks。两个 smoke 作业退出后共享临时目录
均为空。

stderr 中只有 Apptainer 用户命名空间导致 CMA 不可用、OpenMPI 回退其他
节点内传输机制的性能提示，没有正确性错误。基线跨节点路径使用 TCP；
vendor fabric 加速未在本次验证范围内。

## 9. 失败门槛与闭环

首轮 SIF 构建 `60088222` 和完整编译 `60088635` 成功，但单节点 smoke
`60088685` 在 `MPI_Init` 失败。根因不是宿主/容器 MPI 版本不一致，而是
OpenMPI 虽然具有 `ess:pmi`，构建时却显式禁用了 Slurm 环境识别组件。

修复是在同一个 Dockerfile、同一个 `runtime` target 和同一个 SIF 中同时
启用 `--with-pmi=/opt/pmi2` 与 `--with-slurm`，并增加
`MCA schizo: slurm` 构建断言。旧 Docker archive、SIF 和 runner 以带
`noslurm-60088685` 的名称保留为失败证据；没有增加 pilogin 专属镜像，也
没有回退到宿主 `mpirun`。

逐项原始证据索引见 `validation/pilogin-20260730.md`。
