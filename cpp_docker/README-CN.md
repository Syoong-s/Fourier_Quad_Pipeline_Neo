# Fourier_Quad C++ 容器

本目录构建一个供 `cpp_Standard` 或 `cpp_Lite` 使用的 x86_64 Linux 工具链镜像。
流水线源码、配置、星表、观测数据和输出不写入镜像，始终由宿主 bind 挂载。

> English: [README.md](README.md)

## 运行环境

镜像基于 Rocky Linux 8.10，包含 G++ 12.3.0、OpenMPI 4.1.8（PMI2）、
CFITSIO 4.6.4、FFTW 3.3.11、Eigen 3.4.0、LAPACK 3.11.0 与
OpenBLAS 0.3.33。

通用 HPC 基线要求 x86_64、Slurm `pmi2`、Apptainer/Singularity、所有节点同路径
可见的共享文件系统和可路由 TCP。只提供 PMIx、ARM、其他调度器或 vendor fabric
加速的站点需要另行验证。

## 构建与验证

```bash
docker build --platform linux/amd64 --target runtime \
  --build-arg BUILD_JOBS=4 \
  -t cpppipeline-dev:gxx12.3-openmpi4.1.8-pmi2 .
bash scripts/verify-image.sh cpppipeline-dev:gxx12.3-openmpi4.1.8-pmi2
```

验证脚本检查版本、科学库、双 MPI rank、PMI 支持，并确认镜像不含流水线源码。

## 本地使用

```bash
cp .env.example .env
# 设置 CPP_SOURCE_HOST，以及本次启用阶段需要的宿主路径。
docker compose run --rm FourierQuad-CPP
```

容器内执行：

```bash
make -C /workspace/src_pipe -j4
cp /workspace/src_pipe/pipeline.example.ini \
   /workspace/src_pipe/pipeline.ini
# pipeline.ini 必须写容器路径，而不是宿主路径。
mpirun -np 4 /workspace/src_pipe/Fourier_Quad_Pipe \
  --config /workspace/src_pipe/pipeline.ini
```

核心挂载为源码、测天/源星表、平场和可写处理目录。Science/DQ 归档以及
extcat/rearr/曝光表/FD 挂载按阶段选用；Docker Compose 中仅在需要时叠加
`compose.optional.yaml`。

当前程序可在 `pipeline.ini` 设置测天/源星表路径和其他运行期参数。仍位于
`config/*.hpp` 的固定数值参数需要重新编译。

运行 `process_astrocat` 时，应通过合适的只读 bind 暴露原始 Gaia 目录，并在
`[astrocat].input_directory` 或 `--astrocat-input` 中使用其容器路径。星表 bind 是
只读的，因此须把 `[astrocat].output_directory` 或 `--astrocat-output` 指到可写位置，
通常放在处理数据 bind 下。该生产端输出与 `[lensing].astrometry_cat` 相互独立；
后续消费作业须另行让后者指向生成目录，并设置
`[lensing].astrometry_cat_type = 2`。

## Slurm

同一镜像可转换为一个 SIF。复制 `runner/cpppipeline.env.example`，再按
[runner 中文指南](runner/README-CN.md) 操作。支持的启动链为：

```text
srun --mpi=pmi2 -> run-apptainer.sh -> apptainer exec --cleanenv -> Fourier_Quad_Pipe
```

依赖来源与许可见 [SOURCES.md](SOURCES.md) 和
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)。
