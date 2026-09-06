# Fourier_Quad C++ 容器

本目录构建供 `cpp_Standard` 或 `cpp_Lite` 使用的 x86_64 Linux 工具链镜像。
源码、配置、星表、观测数据和输出位于镜像外，通过宿主 bind 挂载。

> English: [README.md](README.md)

## 运行环境

Rocky Linux 8.10 镜像包含 G++ 12.3.0、OpenMPI 4.1.8（PMI2）、
CFITSIO 4.6.4、FFTW 3.3.11、Eigen 3.4.0、LAPACK 3.11.0 和
OpenBLAS 0.3.33。

通用 HPC 基线要求 x86_64、Slurm `pmi2`、Apptainer/Singularity、共享文件系统和
可路由 TCP。其他架构、只提供 PMIx 的站点、其他调度器或 vendor fabric 需另行验证。

## 构建与验证

### 拉取 GHCR 镜像

```bash
docker pull ghcr.io/syoong-s/fourier_quad_pipeline_neo:latest
```

### 下载源码并构建

从 [GitHub Releases](https://github.com/Syoong-s/Fourier_Quad_Pipeline_Neo/releases)
下载 `cpp-docker-<tag>.zip`，解压后执行：

```bash
docker build --platform linux/amd64 --target runtime \
  --build-arg BUILD_JOBS=4 \
  -t cpppipeline-dev:gxx12.3-openmpi4.1.8-pmi2 .
bash scripts/verify-image.sh cpppipeline-dev:gxx12.3-openmpi4.1.8-pmi2
```

验证脚本检查组件版本、科学库、双 MPI rank、PMI 支持，并确认镜像不含 Pipeline 源码。

## 本地使用

```bash
cp .env.example .env
# 设置 CPP_SOURCE_HOST 和启用阶段需要的宿主路径。
docker compose run --rm FourierQuad-CPP
```

### 常改 `.env` 参数

复制 `.env.example` 后按宿主目录修改下列值。`*_HOST` 是宿主路径；对应的
`*_CONTAINER` 是容器内 `pipeline.ini` 或 CLI 使用的绝对路径。

| 参数 | 通常如何修改 | 约束 |
|---|---|---|
| `IMAGE_NAME` | 设为准备运行或本地构建的镜像 tag。 | 必须与已拉取镜像或 `docker build -t` 一致。 |
| `BUILD_JOBS` | 设置镜像构建并行度。 | 按可用 CPU 和内存调整。 |
| `HOST_UID`、`HOST_GID` | 设置当前宿主用户 UID/GID。 | 输出需要由宿主用户直接读写时设置。 |
| `CPP_SOURCE_HOST` | 指向 `cpp_Standard` 或 `cpp_Lite` 源码目录。 | 以读写方式挂载到 `/workspace/src_pipe`，以保留配置和编译产物。 |
| `SCIENCE_ROOT_*`、`DQ_ROOT_*` | 指向 Science/DQ 归档。 | 需要时通过 `compose.optional.yaml` 启用；Lite 主流程始终需要 DQ。 |
| `ASTROMETRY_CAT_*`、`SOURCE_CAT_*`、`FLAT_PATH_*` | 指向 Gaia、External source catalog 和平场目录。 | 容器路径必须与有效 INI/CLI 值一致；flat 仅相应 Standard 分支需要。 |
| `PROCESS_DATA_*` | 指向可写处理目录。 | 保存曝光表、中间文件和结果；容器内默认 `/data/DataProcess`。 |
| `EXTCAT_INPUT_*`、`REARR_OUTPUT_*`、`EXPOLIST_DIR_*`、`FD_OUTPUT_*` | 设置可选的独立阶段路径。 | 使用时叠加 `compose.optional.yaml`；否则优先放在处理目录下。 |

容器内执行：

```bash
make -C /workspace/src_pipe -j4
cp /workspace/src_pipe/pipeline.example.ini \
   /workspace/src_pipe/pipeline.ini
# 在 pipeline.ini 中填写容器路径。
mpirun -np 4 /workspace/src_pipe/Fourier_Quad_Pipe \
  --config /workspace/src_pipe/pipeline.ini
```

核心 bind 为源码、测天/源星表、平场和可写处理目录。Science/DQ 归档及独立阶段路径
按所选阶段挂载。运行时路径写入 `pipeline.ini` 或受支持的 CLI；固定
`config/*.hpp` 值修改后需要重新编译。

运行 `process_astrocat` 时，应只读 bind 原始 Gaia 目录，并把分片写到可写路径，
通常放在处理数据 bind 下。该生产端输出与 `[lensing].astrometry_cat` 相互独立；
消费路径须另行配置。Standard 还须设置 `[lensing].astrometry_cat_type = 2`；Lite 已固定
使用一度分片布局。

## Slurm

将同一镜像转换为一个 SIF，然后按 [runner 中文指南](runner/README-CN.md) 操作。
支持的启动链为：

```text
srun --mpi=pmi2 -> run-apptainer.sh -> apptainer exec --cleanenv -> Fourier_Quad_Pipe
```

依赖来源与许可见 [SOURCES.md](SOURCES.md) 和
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)。
