# Fourier_Quad C++ Slurm/Apptainer Runner

本 runner 在提供 `pmi2` 的 x86_64 Slurm 集群上启动 SIF 内链接的流水线。编译器、
OpenMPI 与科学库由 SIF 提供，不要求宿主编译器/OpenMPI ABI 与应用一致。

> English: [README.md](README.md)

## 前提

- `srun --mpi=list` 包含 `pmi2`；
- 计算节点可用 Apptainer 或 Singularity；
- runner、SIF、源码、星表、数据、cache、tmp 位于所有节点同路径可见的共享存储；
- 分配节点之间可通过 TCP 通信。

先做只读检查：

```bash
bash inspect-cluster-mpi.sh
```

## 配置

```bash
cp cpppipeline.env.example cpppipeline.env
```

设置 SIF 或 OCI/archive 来源、`CPP_SOURCE_HOST`、`PROCESS_DATA_HOST`、星表/平场
bind 和 Apptainer cache/tmp。Science、DQ、extcat 输入、rearr 输出、曝光表、FD 输出
只在所选阶段需要时设置。

`cpppipeline.env` 会被 Bash `source`。`HPC_MODULES`、`HPC_EXTRA_BINDS`、
`HPC_PASSTHROUGH_ENV`、`HPC_CONTAINER_ENV`、`SRUN_ARGS` 必须保持数组，并保留：

```text
MPI_LAUNCH_MODE=srun
SLURM_MPI_TYPE=pmi2
```

模块可以提供 Slurm/Apptainer，但不得把宿主 MPI 库注入应用环境。

## 获得 SIF

从已审核 Docker archive 构建：

```bash
sbatch build-sif.slurm
```

从固定 digest 的 registry 镜像拉取：

```bash
bash pull-sif.sh
```

两种方式都拒绝覆盖已有 SIF，并生成 SHA256 sidecar。生产运行前将其值写入
`CPP_SIF_SHA256_EXPECTED`。

## 验证与运行

```bash
bash run-apptainer.sh --check
sbatch compile-pipeline.slurm
sbatch --nodes=1 --ntasks=2 --ntasks-per-node=2 mpi-smoke-test.slurm
sbatch mpi-smoke-test.slurm
```

全部通过后再提交真实流水线：

```bash
sbatch cpppipeline.slurm \
  --config /workspace/src_pipe/pipeline.ini
```

脚本名后的参数会原样转发给 `Fourier_Quad_Pipe`；没有额外参数时，runner 将
`CPP_EXPO_LIST_CONTAINER` 作为兼容位置参数传入。

源码 bind 可写，编译产物保留在其中；不要同时编译同一源码副本。星表/平场只读，
处理与输出目录可写。`.slurm` 中的 partition、account、节点、task、CPU、内存、时间和
日志都是模板，应按站点整体调整，但不要改变 `srun --mpi=pmi2` 启动边界。

运行 `process_astrocat` 时，在 `HPC_EXTRA_BINDS` 中把原始 Gaia 目录设为只读 bind，
再把生成分片写到可写的 `PROCESS_DATA_CONTAINER` 下或另一个明确的可写 bind。
在 `[astrocat]` 或对应 CLI 选项中使用这些容器路径。生产端输出与
`[lensing].astrometry_cat` 相互独立；消费路径和
`[lensing].astrometry_cat_type = 2` 须另行配置。
