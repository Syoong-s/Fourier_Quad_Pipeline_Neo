# cpp_Standard 通用 Slurm/Apptainer runner

本 runner 使用一个预编译 SIF，面向所有满足运行契约的 x86_64 Slurm 集群。
SIF 自带 G++ 12.3.0、OpenMPI 4.1.8、PMI2 client 和科学库；宿主
GCC/OpenMPI 不进入应用 ABI。

启动链只有一条：

`srun --mpi=pmi2` → `run-apptainer.sh` → `apptainer exec --cleanenv` →
SIF 内链接的 `Fourier_Quad_Pipe`。

## 1. 站点前提

- x86_64 Linux 计算节点；
- `srun --mpi=list` 包含 `pmi2`；
- 计算节点可使用 Apptainer 或 Singularity；
- SIF、runner、源码、数据和 scratch 位于所有节点同路径可见的共享文件系统；
- 分配节点之间有可路由 TCP。

基线镜像没有绑定某个站点的 UCX、OFI 或 RDMA provider。它验证可运行性，
不把 vendor fabric 性能当作通用能力。

## 2. 推荐目录

```text
/shared/project/cpppipeline/
├── code/                  # 独立、可写的 cpp_Standard 源码副本
├── runner/                # 本目录完整副本
├── images/                # Docker archive、SIF、SHA256 sidecar
├── apptainer-cache/
├── apptainer-tmp/
├── scratch/               # 所有计算节点可写
└── data/
    ├── Science/               # 原始压缩 Science 归档，只读
    ├── DQMask/                # 原始压缩 DQMask 归档，只读
    ├── AstroDir/
    ├── ExtSrcDir/
    ├── FlatDir/
    └── DataProcess/
```

不要让两个作业同时编译同一份 `code/`。

## 3. 审计和配置

先运行只读站点审计：

```bash
bash inspect-cluster-mpi.sh
```

确认输出包含 `pmi2`，然后复制配置：

```bash
cp cpppipeline.env.example cpppipeline.env
```

编辑所有宿主绝对路径。以下变量是 Bash indexed array，不能改成普通字符串：

- `HPC_MODULES`
- `HPC_EXTRA_BINDS`
- `HPC_PASSTHROUGH_ENV`
- `HPC_CONTAINER_ENV`
- `SRUN_ARGS`

`HPC_MODULES` 可以初始化 Apptainer/Slurm，但不应加载宿主 MPI。
`run-apptainer.sh` 使用 `--cleanenv` 清除宿主环境污染，再转发 Slurm 创建的
全部 `SLURM_*`、`PMI_*`、`PMI2_*`，以及 allowlist 中的额外变量。必须保持
`MPI_LAUNCH_MODE=srun`、`SLURM_MPI_TYPE=pmi2` 和
`HPC_SCRUB_MPI_ENV=1`。

`SCIENCE_ROOT_HOST` 和 `DQ_ROOT_HOST` 是可选的只读归档挂载源。仅在运行
process_init 或 process_main 时设置。对应容器路径为
`SCIENCE_ROOT_CONTAINER` 和 `DQ_ROOT_CONTAINER`。四个额外的可选挂载
（`EXTCAT_INPUT`、`REARR_OUTPUT`、`EXPOLIST_DIR`、`FD_OUTPUT`）遵循相同的
条件模式：设置 `*_HOST` 即绑定，留空则跳过。生成结果继续写入可写的
`PROCESS_DATA_CONTAINER`。

三个 catalogue/flat 容器路径必须与
`/workspace/src_pipe/include/process_main/LensingConfig.hpp` 编译期常量一致。

## 4. 获得 SIF

从已审核 Docker archive 构建时，配置 `CPP_DOCKER_ARCHIVE`、`CPP_SIF`、
`APPTAINER_CACHE_DIR` 和 `APPTAINER_TMP_DIR`，提交：

```bash
sbatch build-sif.slurm
```

从 registry 拉取时，把 `OCI_IMAGE_URI` 固定为 digest，执行：

```bash
bash pull-sif.sh
```

两个入口都先创建同目录临时文件，成功后原子改名，同时生成
`${CPP_SIF}.sha256`，并拒绝覆盖已有文件。正式运行前将 sidecar 第一列填入
`CPP_SIF_SHA256_EXPECTED`。

## 5. 单实例检查

```bash
bash run-apptainer.sh --check
```

该检查在宿主侧验证可选 SIF hash，在容器侧验证唯一 image ID、G++、
OpenMPI、`ess:pmi`、`libpmi2`、CFITSIO、FFTW、Eigen 和所有基础 bind。

## 6. 编译完整 pipeline

```bash
sbatch compile-pipeline.slurm
```

此作业只有一个 task，在 SIF 中按 `CPP_MAKE_CLEAN` 执行清理，再使用
`CPP_BUILD_JOBS` 并行编译 22 个 translation unit，最终确认
`Fourier_Quad_Pipe` 可执行。`CPP_BUILD_JOBS` 不得超过申请的
`SLURM_CPUS_PER_TASK`。

## 7. PMI2 烟雾测试

先运行单节点双 rank：

```bash
sbatch \
    --nodes=1 \
    --ntasks=2 \
    --ntasks-per-node=2 \
    mpi-smoke-test.slurm
```

再运行模板默认的双节点四 rank：

```bash
sbatch mpi-smoke-test.slurm
```

作业在 `HPC_SHARED_SCRATCH_HOST` 创建唯一临时目录，编译 MPI identity 和
科学库测试，再由 `srun --mpi=pmi2` 启动每个 SIF rank。成功日志应包含
全部 rank、预期 hostname 和科学栈成功消息。退出 trap 只清理该作业自己
创建的目录。

## 8. 运行真实 pipeline

确认 exposure list 和编译期科学路径后提交：

```bash
sbatch cpppipeline.slurm
```

无额外参数时，程序读取
`${PROCESS_DATA_CONTAINER}/expo_list.list`。脚本名后的参数会传给
`Fourier_Quad_Pipe`。真实数据运行与无数据 smoke 分开，runner 不会自行
修改科学路径或 exposure 内容。

串联运行时，在脚本名后传入 `--run-init true --run-main true`、
`--science-root /data/archive/science`、`--dq-root /data/archive/dqmask`、
`--output-root /data/DataProcess`，并可重复传入例如
`--dataset g2013:c4d_13 --dataset g2014:c4d_14`。多个 `--dataset` 会顺序
执行；多个 `--contains` 按 OR 规则匹配，runner 会原样转发这些参数。

## 9. Slurm 模板与日志

三个 `.slurm` 文件的 partition、account、节点、task、CPU、内存和时间都
是模板，应按站点政策成组覆盖。若站点要求集中日志，每次提交都通过
`sbatch --output=绝对路径 --error=绝对路径` 指定位置。

文件职责：

- `cpppipeline.env.example`：唯一通用配置模板；
- `build-sif.slurm` / `pull-sif.sh`：获得唯一 SIF；
- `run-apptainer.sh`：统一 bind、环境隔离和容器入口；
- `compile-pipeline.slurm`：单 task 完整编译；
- `mpi-smoke-test.slurm`：单/多节点 PMI2 和科学库 smoke；
- `cpppipeline.slurm`：只启动已编译的真实 pipeline。
