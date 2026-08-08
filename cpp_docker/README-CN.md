# cpp_Standard 单 SIF 可移植 HPC 环境

本目录参考 `f77_docker` 已验证的制品边界：镜像只保存工具链与依赖，
`cpp_Standard` 源码、压缩 Science/DQ 归档、星表、平场、处理数据和输出始终
在镜像外，通过 bind 挂载。区别在于这里不再为本地和某个集群分别构建镜像，
而是只构建一个 x86_64 Linux OCI 镜像，并转换成一个预编译 SIF。

## 1. 唯一运行时契约

| 组件 | 版本或接口 |
| --- | --- |
| G++ | 12.3.0 |
| OpenMPI | 4.1.8 |
| Slurm 客户端接口 | Slurm 25.11.2 的 PMI2 |
| OpenMPI Slurm 集成 | 启用 direct-launch 环境识别组件 |
| CFITSIO | 4.6.4 |
| FFTW | 3.3.11 |
| Eigen | 3.4.0 |
| LAPACK / BLAS | 3.11.0 |
| OpenBLAS | 0.3.33 |

项目只有：

- 一个 Docker target：`runtime`；
- 一个 Pixi environment：`default`；
- 一个 image ID：`gxx12.3-openmpi4.1.8-pmi2`；
- 一个供各站点复用的 SIF。

集群执行时不加载宿主 GCC/OpenMPI，不调用宿主 `mpirun`，也不要求宿主 MPI
与镜像 MPI 版本一致。Slurm 通过 `srun --mpi=pmi2` 为每个 task 建立 PMI2
状态，每个 task 再进入同一个 SIF；应用只链接 SIF 中的 OpenMPI 和科学库。

“可移植 HPC”在本项目中的精确定义是：x86_64 Linux、Slurm 提供 PMI2
plugin、计算节点可运行 Apptainer 或 Singularity，并且相关路径位于共享
文件系统。它不自动覆盖 ARM、非 Slurm 调度器或只提供 PMIx 的站点。
基线通信为节点内共享内存和节点间可路由 TCP；UCX、OFI、RDMA 等加速路径
必须按站点另行验证，但不应因此派生站点专属基础镜像。

## 2. 构建和本地验证

进入本目录执行：

```bash
docker build \
    --platform linux/amd64 \
    --target runtime \
    --build-arg BUILD_JOBS=4 \
    -t cpppipeline-dev:gxx12.3-openmpi4.1.8-pmi2 \
    .

bash scripts/verify-image.sh \
    cpppipeline-dev:gxx12.3-openmpi4.1.8-pmi2

bash scripts/check-public-repo.sh
```

验证脚本检查 G++、OpenMPI、`ess:pmi`、`libpmi2`、科学库和镜像标签，
编译并运行两个 MPI rank 的 identity 与科学库测试，并确认最终镜像没有
pipeline 源码。

## 3. 本地挂载与完整编译

复制配置并填写宿主绝对路径：

```bash
cp .env.example .env
docker compose run --rm FourierQuad-CPP
```

容器内执行：

```bash
make -C /workspace/src_pipe clean
make -C /workspace/src_pipe -j4
```

Makefile 使用 `mpicxx`、C++17、CFITSIO、FFTW、Eigen、LAPACK 和 BLAS。
可选的 `STACK_PREFIX` 用于非容器环境。编译产物留在宿主挂载的源码目录。

必需容器路径为：

- 源码：`/workspace/src_pipe`
- 测天星表：`/data/catalogs/AstroDir`
- 源星表：`/data/catalogs/ExtSrcDir`
- 平场：`/data/calib/FlatDir`
- 处理数据：`/data/DataProcess`

六个可选挂载（Science 归档、DQMask、extcat 输入、rearr 输出、expo list、
FD 输出）可通过 `compose.optional.yaml` 引入：

    docker compose -f compose.yaml -f compose.optional.yaml up

HPC（Apptainer）中，在 `cpppipeline.env` 设置对应的 `*_HOST` 变量即可绑定；
留空则跳过。

三个 catalogue/flat 路径必须与
`cpp_Standard/config/LensingConfig.hpp` 的编译期科学配置一致。
初始化时，将 `SCIENCE_ROOT_CONTAINER`、`DQ_ROOT_CONTAINER` 分别传给
`--science-root`、`--dq-root`（仅当这些挂载激活时），并将可写的
`PROCESS_DATA_CONTAINER` 传给 `--output-root`。容器挂载接口无需为批处理
增加变量：重复传入 `--dataset TARGET:PREFIX` 即可批量运行，重复传入
`--contains TOKEN` 即按 OR 规则匹配。

## 4. 生成唯一 SIF

可以把已经验证的镜像保存为 Docker archive：

```bash
docker save \
    --output cpppipeline-gxx12.3-openmpi4.1.8-pmi2.docker.tar \
    cpppipeline-dev:gxx12.3-openmpi4.1.8-pmi2
```

在集群共享目录中复制 `runner/cpppipeline.env.example` 为
`runner/cpppipeline.env`，配置 `CPP_DOCKER_ARCHIVE`、`CPP_SIF`、
Apptainer cache/tmp 和其他共享路径，然后在计算节点提交
`runner/build-sif.slurm`。

如果站点可以访问 OCI registry，也可把 `OCI_IMAGE_URI` 固定为 digest 后
运行 `runner/pull-sif.sh`。两个入口都先写临时文件，成功后原子改名，
生成 `.sha256` sidecar，并拒绝覆盖已有 SIF。

## 5. 集群验证顺序

集群上按以下顺序执行：

1. `run-apptainer.sh --check`：校验 SIF hash、精确版本、PMI 组件和 bind；
2. `compile-pipeline.slurm`：单 task 编译完整 pipeline；
3. `mpi-smoke-test.slurm`：先单节点双 rank，再多节点；
4. `cpppipeline.slurm`：仅在科学路径和 exposure list 复核后运行真实数据。

runner 使用 `--cleanenv` 阻断宿主 MPI 环境污染，再显式转发 Slurm 生成的
`SLURM_*`、`PMI_*`、`PMI2_*`，以及 env 文件中的显式 allowlist。站点可通过
`SITE_ENV_SCRIPT` 和 `HPC_MODULES` 初始化 Apptainer，但不得把宿主 MPI
加入应用 ABI 路径。

OpenMPI 构建同时启用 PMI 和 Slurm direct-launch 组件。后者只负责识别
Slurm 环境与直接启动状态，不把应用链接到宿主 `libslurm`。

完整操作见 [runner/README-CN.md](runner/README-CN.md)。架构决策与验证记录
见 [DEVELOPMENT-SUMMARY-CN.md](DEVELOPMENT-SUMMARY-CN.md)，第三方来源与
许可边界见 [SOURCES.md](SOURCES.md) 和
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)。
