# Fourier_Quad C++ Slurm/Apptainer Runner

本 runner 在提供 `pmi2` 的 x86_64 Slurm 集群上启动 SIF 内链接的 Pipeline。
编译器、OpenMPI 与科学库由 SIF 提供，不要求宿主编译器/OpenMPI ABI 与应用一致。

> English: [README.md](README.md)

## 前提与配置

确认 `srun --mpi=list` 包含 `pmi2`、计算节点可用 Apptainer/Singularity、分配
节点间可通过 TCP 通信，且所有路径在各节点同位置可见。然后执行：

```bash
bash inspect-cluster-mpi.sh
cp cpppipeline.env.example cpppipeline.env
```

### 常改 runner 参数

`cpppipeline.env` 中的宿主路径必须在所有分配节点共享；成对的容器路径必须与
`pipeline.ini` 或对应 CLI 使用的路径一致。

| 参数或位置 | 通常如何修改 | 约束 |
|---|---|---|
| `OCI_IMAGE_URI` / `CPP_DOCKER_ARCHIVE` | 按 OCI 拉取或本地 Docker archive 转换二选一设置。 | URI 应固定 digest；archive 必须是构建作业可见的 x86_64 镜像。 |
| `CPP_SIF`、`CPP_SIF_SHA256_EXPECTED` | 设为共享 SIF 路径和可选预期 SHA256。 | 所有节点使用同一路径；获得镜像后再填写校验值。 |
| `CPP_SOURCE_HOST/CONTAINER` | 指向 `cpp_Standard` 或 `cpp_Lite`，通常挂载到 `/workspace/src_pipe`。 | 宿主源码需可写以保存配置和编译产物。 |
| `SCIENCE_ROOT_*`、`DQ_ROOT_*` | 指向初始化或主流程使用的 Science/DQ 归档。 | 只在需要时 bind；Lite 主流程需要逐 CCD DQ。 |
| `ASTROMETRY_CAT_*`、`SOURCE_CAT_*`、`FLAT_PATH_*` | 指向 Gaia、External source catalog 和平场目录。 | 容器路径须与有效 INI/CLI 值一致；flat 仅 Standard 使用。 |
| `PROCESS_DATA_*`、`CPP_EXPO_LIST_CONTAINER` | 设为共享可写处理目录和默认曝光表路径。 | 曝光表和配置必须位于已 bind 的路径下。 |
| `EXTCAT_INPUT_*`、`REARR_OUTPUT_*`、`EXPOLIST_DIR_*`、`FD_OUTPUT_*` | 设置可选的独立阶段目录。 | 未设置时不 bind；INI/CLI 必须使用相应容器路径。 |
| `HPC_SHARED_SCRATCH_HOST`、`APPTAINER_CACHE_DIR`、`APPTAINER_TMP_DIR` | 设为共享可写 scratch/cache/tmp。 | 应有足够空间且计算节点可访问。 |
| `APPTAINER_BIN`、`HPC_MODULES`、`SITE_ENV_SCRIPT` | 按站点运行时和模块环境设置。 | 不得注入宿主 MPI 库。 |
| `MPI_LAUNCH_MODE=srun`、`SLURM_MPI_TYPE=pmi2`、`SRUN_ARGS=()` | 保留启动方式，必要时通过 `SRUN_ARGS` 增加站点参数。 | 本 runner 要求 `srun` + `pmi2`；保留数组语法。 |
| `HPC_EXTRA_BINDS=()`、`HPC_PASSTHROUGH_ENV=(OMP_NUM_THREADS)`、`HPC_CONTAINER_ENV=()` | 只设置必要的额外文件或环境变量。 | 保留 Bash 数组，并最小化传入 `--cleanenv` 的宿主状态。 |
| `CPP_BUILD_JOBS`、`CPP_MAKE_CLEAN`、`CPP_EXECUTABLE` | 调整编译并行度、是否先清理和可执行文件路径。 | 不要让多个作业同时清理或编译同一源码副本。 |
| 镜像/版本预期值 | 更换镜像软件栈时同步修改。 | 必须与 SIF 内实际版本一致。 |
| `cpppipeline.slurm` 中的 `#SBATCH` | 修改站点 partition、account、节点、task、CPU、内存、时间和日志。 | 遵守站点策略，并保持 task 数与 MPI rank 一致。 |

`cpppipeline.env` 是 Bash；保留其中的数组，并保持：

```text
MPI_LAUNCH_MODE=srun
SLURM_MPI_TYPE=pmi2
```

模块可提供 Slurm/Apptainer，但不得把宿主 MPI 库注入应用。

## 获得与验证 SIF

任选一种来源：

```bash
sbatch build-sif.slurm          # 本地已有 Docker archive
bash pull-sif.sh                # 远程 OCI 镜像
```

二者都拒绝覆盖已有 SIF，并生成 SHA256 sidecar。生产运行前将其值写入
`CPP_SIF_SHA256_EXPECTED`，随后执行：

```bash
bash run-apptainer.sh --check
sbatch compile-pipeline.slurm
sbatch --nodes=1 --ntasks=2 --ntasks-per-node=2 mpi-smoke-test.slurm
sbatch mpi-smoke-test.slurm
```

## 运行

在脚本名后传入 INI 路径或 C++ CLI：

```bash
sbatch cpppipeline.slurm \
  --config /workspace/src_pipe/pipeline.ini
```

参数会原样转发给 `Fourier_Quad_Pipe`；没有参数时，runner 将
`CPP_EXPO_LIST_CONTAINER` 作为兼容曝光表位置参数传入。源码 bind 可写以保存编译
产物；不要并发编译同一副本。按站点修改 Slurm 资源，但不要改变
`srun --mpi=pmi2` 启动边界。

运行 `process_astrocat` 时，在 `HPC_EXTRA_BINDS` 中把原始 Gaia 目录设为只读 bind，
并把生成分片写到可写 bind 下。在 `[astrocat]` 或 CLI 中配置生产端路径；独立的
`[lensing].astrometry_cat` 消费路径须另行配置。Standard 还须设置
`[lensing].astrometry_cat_type = 2`；Lite 已固定使用一度分片布局。
