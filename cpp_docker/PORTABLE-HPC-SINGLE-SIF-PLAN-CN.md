# cpp_Standard 单一 SIF 跨 HPC Docker 与 Runner 实施方案

## 1. 结论

`cpp_Standard` 只构建和发布一个预编译 SIF，不再区分 `local`、
`pilogin` 或其他集群镜像。

镜像内固定编译器、MPI 和科学计算依赖；集群只负责资源分配、启动
容器 rank、挂载代码与数据。Slurm 集群的标准启动路径为：

```text
sbatch
  └─ srun --mpi=pmi2
       ├─ apptainer exec <同一个 SIF> <rank 0>
       ├─ apptainer exec <同一个 SIF> <rank 1>
       └─ ...
```

应用不使用宿主 `mpirun`，不链接宿主 OpenMPI，也不要求宿主 OpenMPI
与镜像内 OpenMPI 同版本。pilogin 只是该通用制品的一个验证环境，
不是镜像 target、tag、profile 或 runner 分支。

## 2. 目标与可移植性边界

### 2.1 必须达到

- 同一份 x86_64 Linux SIF 原样复制到不同 HPC，不重新编译镜像。
- `cpp_Standard` 始终使用 SIF 内的 `mpicxx` 编译并链接。
- 首个正式运行契约是：
  - Slurm 提供 `srun`；
  - `srun --mpi=list` 包含 `pmi2`；
  - 计算节点可以执行 Apptainer 或 Singularity；
  - SIF、源码和运行数据位于所有分配节点可见的共享文件系统。
- runner 只适配资源参数、路径、容器运行时模块和站点初始化。
- 单节点、多节点使用同一套脚本和同一个 SIF。
- 构建源、版本、校验和、镜像 digest、SIF SHA-256 和集群测试结果均可追溯。

### 2.2 不作不真实承诺

“任意 HPC”不能理解为不检查架构、调度器或进程管理接口便无条件运行。
第一版保证的是“满足上述 Slurm/PMI2 契约的 x86_64 HPC”。以下场景不
通过新增站点镜像解决：

- 集群没有 PMI2：先失败并给出清晰诊断；PMIx 只能作为经过独立验证的
  显式 runner 模式，不能静默回退。
- PBS、LSF 等非 Slurm 调度器：以后增加 launcher adapter，但仍复用
  同一个 SIF。
- ARM 或其他 CPU 架构：需要同一配方的对应架构制品，不可能让一个
  x86_64 SIF 直接运行。
- InfiniBand、UCX、OFI 等高性能网络：基础功能以共享内存和 TCP 为准；
  原生互连加速需要逐集群验证，但不能反向污染通用镜像设计。

## 3. 从 f77 runner 继承的结构

f77 的已验证实现提供了四个关键原则：

| f77 已验证做法 | cpp_Standard 对应设计 |
| --- | --- |
| 一个固定的预编译 SIF | 一个固定的 cpp SIF |
| 应用链接容器内 MPI | 应用链接容器内 OpenMPI 4.1.8 |
| `srun --mpi=pmi2` 直接启动每个容器 rank | 使用相同的 Slurm/PMI2 启动拓扑 |
| runner 保留 Slurm/PMI 环境并隔离宿主 MPI 污染 | cpp wrapper 执行同样的选择性环境清理 |
| 编译只进行一次，MPI 阶段再启动多个 rank | `compile-pipeline.slurm` 与运行 job 分离 |
| pilogin 是验证记录 | pilogin 只产生审计记录，不产生专属文件 |

因此本方案不采用 host/container OpenMPI hybrid，也不要求加载
`gcc/12.3.0 openmpi/4.1.6-gcc-12.3.0` 之类的宿主编译器和 MPI 模块。

## 4. 唯一镜像基线

### 4.1 制品标识

建议使用与站点无关的名称：

```text
OCI tag: cpppipeline:gxx12.3-openmpi4.1.8-pmi2
SIF:     cpppipeline-gxx12.3-openmpi4.1.8-pmi2.sif
```

发布时必须同时记录不可变 OCI digest 和 SIF SHA-256。runner 的生产
配置以 digest/校验和识别制品，不以可变 tag 作为唯一依据。

### 4.2 固定软件版本

| 组件 | 计划版本 | 说明 |
| --- | --- | --- |
| 基础用户态 | Rocky Linux 8.10 | 采用较保守的 glibc 基线；实现时固定镜像 digest |
| GCC/G++ | 12.3.0 | 与现有 cpp pilogin 配置一致，编译器随 SIF 分发 |
| OpenMPI | 4.1.8 | OpenMPI 4.1 稳定分支版本，启用外部 PMI2 |
| CFITSIO | 4.6.4 | 复用当前锁定版本 |
| FFTW | 3.3.11 | 复用当前锁定版本 |
| Eigen | 3.4.0 | 复用当前锁定版本 |
| BLAS/LAPACK | 3.11.0 | OpenBLAS 实现 |
| OpenBLAS | 0.3.33 | 复用当前锁定版本 |

OpenMPI 4.1.8 源码使用官方归档：

```text
https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.8.tar.bz2
SHA-256:
466f68e3132a1dc02710cc2011fafced8336d98359fa2dae4dddcfd5719f12a9
```

### 4.3 PMI2 客户端

OpenMPI 4.1 自身不应在运行时寻找某个站点路径下的 `libpmi2.so`。镜像
构建过程必须：

1. 从带版本和 SHA-256 的官方 Slurm 源码构建 PMI2 客户端头文件及库；
2. 只把 PMI2 客户端运行所需内容安装到 `/opt/pmi2`；
3. 不把 `slurmctld`、`slurmd`、`srun` 或站点配置打入镜像；
4. 使用该固定前缀编译 OpenMPI：

```text
--with-pmi=/opt/pmi2
--without-slurm
--disable-mpi-fortran
--disable-oshmem
--without-cuda
--without-libfabric
--without-psm
--without-psm2
--without-ucx
--without-verbs
```

`--without-slurm` 的含义是镜像不携带或调用 Slurm launcher；直接任务
启动仍由宿主 `srun --mpi=pmi2` 完成。OpenMPI 通过镜像内 PMI2 客户端
连接 Slurm 为当前 job step 提供的 PMI2 服务。

PMI2 客户端的最终 Slurm 源码版本不能仅凭 pilogin 版本决定。实施时先
选择一个仍维护的 Slurm 稳定版本并固定校验和，然后执行第 9 节的跨版本
兼容性门禁。若门禁失败，不允许增加 pilogin 镜像作为补丁；应改选兼容
性更好的 PMI2 客户端版本，或把唯一标准镜像整体切换到 MPICH 路线。

### 4.4 Docker 多阶段结构

保留“多阶段构建”，但最终只有一个 runtime target：

1. `stack-builder`
   - Pixi 只保留一个隐式 `default` environment；
   - 安装 G++ 12.3.0、科学库和构建工具；
   - 不安装 Conda OpenMPI，避免第二套 MPI 混入。
2. `pmi2-builder`
   - 构建固定版本的 Slurm PMI2 client 到 `/opt/pmi2`。
3. `openmpi-builder`
   - 使用同一 G++/GCC 工具链构建 OpenMPI 4.1.8；
   - `--with-pmi=/opt/pmi2`；
   - 安装到 `/opt/openmpi-4.1.8`。
4. `runtime`
   - 复制编译器、科学栈、PMI2 client 和 OpenMPI；
   - 保留 `make`、`pkg-config`、`binutils`，使 HPC 上能编译源码；
   - 不包含 Slurm 命令、SSH daemon、站点模块或站点路径。

最终镜像的 `PATH`、`CPATH`、`LIBRARY_PATH`、`LD_LIBRARY_PATH` 和
`PKG_CONFIG_PATH` 全部指向镜像内固定前缀。构建参数禁止
`-march=native`，采用通用 x86-64 指令基线；OpenBLAS 使用运行时 CPU
分派，避免制品绑定构建机 CPU。

### 4.5 镜像内自检

Docker build 和 SIF `--check` 至少验证：

- `g++ -dumpfullversion` 等于 `12.3.0`；
- `mpicxx --showme:version` 包含 `Open MPI 4.1.8`；
- `ompi_info` 能看到 PMI 支持；
- `ldd`/`readelf` 显示 MPI、PMI2 和 C++ runtime 均来自镜像内前缀；
- CFITSIO、FFTW、Eigen、BLAS/LAPACK 版本和编译链接测试通过；
- 镜像内没有第二套 `mpicc`/`mpicxx` 优先于标准前缀；
- 编译产物不引用 `/usr/lib*/slurm`、宿主模块目录或 pilogin 路径。

## 5. 单一通用 runner

### 5.1 保留和新增的通用文件

实施后 runner 只保留以下与站点无关的接口：

```text
runner/
├── cpppipeline.env.example
├── site-env.example.sh
├── inspect-cluster-mpi.sh
├── pull-sif.sh
├── build-sif.slurm
├── run-apptainer.sh
├── compile-pipeline.slurm
├── mpi-smoke-test.slurm
├── cpppipeline.slurm
└── tests/
    ├── mpi_identity.cpp
    └── science_stack_smoke.cpp
```

以下文件在通用版本验证成功后退役：

```text
compile-pipeline-pilogin.slurm
cpppipeline-pilogin.slurm
cpppipeline.pilogin.env.example
mpi-smoke-test-pilogin.slurm
```

可以保留带日期的 pilogin 测试报告，但测试报告不是 runner 分支。

### 5.2 环境配置契约

`cpppipeline.env.example` 的 MPI 核心值固定为：

```bash
MPI_LAUNCH_MODE=srun
SLURM_MPI_TYPE=pmi2
HPC_SCRUB_MPI_ENV=1
HPC_MODULES=()

CPP_GXX_VERSION_EXPECTED=12.3.0
CPP_OPENMPI_VERSION_EXPECTED=4.1.8
CPP_IMAGE_ID_EXPECTED=gxx12.3-openmpi4.1.8-pmi2
CPP_SIF_SHA256=
```

站点只填写：

- SIF、源码、目录和数据的共享文件系统绝对路径；
- Apptainer/Singularity 可执行文件或负责提供它的模块；
- 可选的站点初始化脚本；
- 分区、账户、QoS、节点、rank、CPU、内存和时间等 Slurm 参数；
- 必要但与 MPI ABI 无关的额外挂载。

`HPC_MODULES` 默认空。若某个集群通过 module 提供 Apptainer，可以只
加载该容器运行时模块；不得为应用执行加载宿主 GCC/OpenMPI 模块。

### 5.3 运行时环境处理

`run-apptainer.sh` 不使用登录 shell，也不全量删除 job 环境。它应：

- 保留 `SLURM_*` 和 PMI2 所需的 `PMI_*`；
- 保留调度器设置的 CPU/GPU binding 信息；
- 去掉宿主 `PATH`、`LD_LIBRARY_PATH`、`LIBRARY_PATH`、`CPATH` 中的
  MPI/编译器模块污染，并由镜像定义这些变量；
- 去掉宿主 `OPAL_PREFIX`、OpenMPI MCA 和 PMIx override；
- 不删除 PMI2 的 `PMI_RANK`、`PMI_SIZE`、`PMI_FD` 等变量；
- 明确挂载源码、目录、校准数据和处理目录；
- 默认 `--no-home`，不依赖用户 home 中的动态库或配置。

环境清理必须是可测试的 allow/preserve 规则，不能用会同时清掉 PMI2
握手信息的粗暴 `env -i`。

### 5.4 启动方式

批处理脚本只允许以下主路径：

```bash
srun --kill-on-bad-exit=1 \
    --ntasks="${SLURM_NTASKS}" \
    --mpi=pmi2 \
    "${RUNNER}" --env-file "${ENV_FILE}" -- \
    "${CPP_EXECUTABLE}" "${pipeline_arguments[@]}"
```

明确禁止：

- `mpirun apptainer exec ...` 的 host/container hybrid；
- 用宿主 OpenMPI launcher 启动镜像内 OpenMPI 应用；
- bind mount 宿主 MPI 安装覆盖镜像 MPI；
- 根据 hostname 或集群名称选择镜像；
- runner 发现 PMI2 不可用后自动改用未经验证的 PMIx。

### 5.5 编译方式

源码和 build 目录位于共享文件系统。`compile-pipeline.slurm` 在 batch
节点上只执行一次：

1. `run-apptainer.sh --check`；
2. 使用镜像内 `make` 和 `mpicxx` 清理并编译；
3. 检查可执行文件；
4. 使用 `ldd`/`readelf` 确认 MPI 与科学库来源；
5. 写入编译器、OpenMPI、git commit 和 SIF SHA-256 到构建审计记录。

编译阶段不启动多个 MPI rank，避免多个 task 同时修改同一 build tree。

## 6. 集群探测与失败前置

`inspect-cluster-mpi.sh` 应输出并保存：

```text
scontrol --version
srun --mpi=list
command -v apptainer || command -v singularity
apptainer version 或 singularity version
共享目录在当前节点的可读写检查
CPU 架构和基础网络接口
```

runner 提交前判定：

1. 未找到 Apptainer/Singularity：失败；
2. `srun --mpi=list` 未出现 `pmi2`：失败并说明当前通用路径不受支持；
3. SIF SHA-256 不匹配：失败；
4. 源码或数据路径不是共享路径/计算节点不可见：失败；
5. 发现宿主 MPI 路径渗入容器：失败；
6. 单 rank 自检未通过：不提交多节点任务。

探测脚本可以报告可用 PMIx 类型，但不自行决定切换。PMIx 支持必须以
独立配置值和完整多节点 smoke 证据启用。

## 7. 从当前实现迁移

| 当前结构 | 目标结构 |
| --- | --- |
| `runtime-local` + `runtime-pilogin` | 单一 `runtime` |
| Pixi `local` + `pilogin` environment | 单一 `default` environment |
| OpenMPI 5.0.10 和 4.1.6 两套 MPI | OpenMPI 4.1.8 + 固定 PMI2 client |
| `local-openmpi5.0.10` / `pilogin-openmpi4.1.6` tag | `gxx12.3-openmpi4.1.8-pmi2` |
| `CPP_STACK_PROFILE_EXPECTED` | 删除 profile，改用唯一 image ID |
| `MPI_LAUNCHER=mpirun` | `MPI_LAUNCH_MODE=srun` |
| 宿主 OpenMPI MCA 参数 | 删除；基础传输由镜像 OpenMPI 管理 |
| `*-pilogin*` runner | 通用 runner + 单独验证报告 |

迁移次序：

1. 合并 `pixi.toml` 与 lockfile，只保留一个软件栈；
2. 重写 Dockerfile 的 PMI2/OpenMPI 构建阶段和单一 runtime；
3. 更新 `verify-image.sh`、Compose 和 OCI labels；
4. 重构 `run-apptainer.sh` 的环境隔离；
5. 把通用 Slurm job 改为 `srun --mpi=pmi2`；
6. 完成镜像、runner 和多节点验证；
7. 只有在新路径全部通过后，删除 pilogin 专属文件；
8. 更新中英文 README、`SOURCES.md` 和第三方许可证说明。

## 8. 构建与分发

### 8.1 唯一制品生成

1. 从 digest 固定的 Rocky 基础镜像构建 OCI；
2. 运行 Docker 级别的单进程和两 rank 测试；
3. 以不可变 OCI digest 生成一次 canonical SIF；
4. 对 SIF 计算 SHA-256；
5. 发布 OCI digest、SIF、SIF SHA-256、SBOM 和源文件校验和；
6. 生产运行优先直接复制该 canonical SIF 并核对其 SHA-256；
7. 若为了验证可复现性而在集群重建 SIF，也必须从同一 OCI digest
   生成，不能使用站点 target。重建 SIF 可能因 Apptainer 版本或元数据
   而具有不同文件 SHA-256，因此应单独记录其校验和，并核对 OCI digest、
   labels、软件版本和镜像内文件清单，不能谎称与 canonical SIF 字节相同。

网络隔离集群可传输 SIF 或 OCI archive。`pull-sif.sh` 和
`build-sif.slurm` 只是分发方式不同，得到的逻辑软件栈必须相同。

### 8.2 pilogin 复现位置

按要求在 pilogin 的 `~/ysx/Agents` 下建立工作目录，例如：

```text
~/ysx/Agents/cpppipeline/
├── image/
├── runner/
├── source/
├── work/
└── validation/
```

该目录中执行：

1. 保存 `inspect-cluster-mpi.sh` 输出；
2. 按原始要求从固定 OCI digest 构建验证用 SIF；同时再放置 canonical
   SIF 时，对后者核对发布的 SHA-256；
3. 对集群重建的 SIF 记录新 SHA-256，并核对 OCI digest、labels、版本
   和镜像内清单；
4. 用通用 `compile-pipeline.slurm` 编译；
5. 运行单节点和双节点 smoke；
6. 最后才提交最小、无破坏性的 pipeline smoke；
7. 将 job ID、脚本 commit、SIF SHA-256、节点/rank 映射和结果写入
   `validation/pilogin-YYYYMMDD.md`。

目录名或验证报告可以包含 pilogin，镜像、runner 和 job 模板不得包含。

## 9. 验证矩阵与准入门禁

### 9.1 镜像级

- Docker 内 `g++`、`mpicxx`、CFITSIO、FFTW、Eigen、BLAS/LAPACK 自检；
- 2 rank 容器内 OpenMPI smoke；
- SIF 单进程 `--check`；
- SIF 内编译两个测试程序；
- 检查 RPATH、动态库来源和第二套 MPI 冲突。

### 9.2 Slurm/PMI2 级

| 场景 | 最低资源 | 通过条件 |
| --- | --- | --- |
| PMI2 启动 | 1 节点、1 rank | MPI 初始化和 finalize 成功 |
| 单节点 MPI | 1 节点、2 ranks | rank 唯一、collective 正确 |
| 多节点 MPI | 2 节点、4 ranks | 4 个唯一 rank，覆盖 2 个 hostname |
| 科学栈 | 2 节点、4 ranks | CFITSIO/FFTW/Eigen/BLAS 测试和 Allreduce 正确 |
| 异常传播 | 2 节点、4 ranks | 任一 rank 失败后 `--kill-on-bad-exit` 使 job 非零退出 |

PMI2 client 兼容性不能只用 pilogin 一次成功来宣称。正式发布前至少在
两个不同 Slurm 主版本、两个计算节点上执行相同 SIF 的 4-rank 测试；
如果当前只能访问 pilogin，应把支持范围标为“pilogin 已验证，其他
Slurm/PMI2 集群待验证”，而不是写成已证明全覆盖。

### 9.3 cpp_Standard 级

- 在 SIF 内完整编译 `cpp_Standard`；
- 编译命令确定使用 `/opt/openmpi-4.1.8/bin/mpicxx`；
- 可执行文件的 MPI、libstdc++、CFITSIO、FFTW、BLAS/LAPACK 均解析到
  SIF 内路径；
- 使用最小测试数据或只读 fixture 完成初始化和至少一个安全处理步骤；
- 所有 rank 正常进入同一 communicator；
- 输出目录使用独立 smoke 路径，不改写生产数据；
- 无真实数据 fixture 时，只完成编译和 MPI/科学栈 smoke，并在报告中
  明确标注 pipeline 数据级 smoke 尚未执行，不能伪报成功。

## 10. 验收标准

只有同时满足以下条件才算完成：

- 仓库中只有一个 cpp runtime target、一个软件环境和一个正式 image tag；
- 没有 `local`/`pilogin` 镜像 profile；
- 通用 runner 默认且明确使用 `srun --mpi=pmi2`；
- host `mpirun`、host OpenMPI 和 host GCC 不在应用 ABI 路径中；
- canonical SIF 具有固定 SHA-256，并能在 pilogin
  `~/ysx/Agents` 中通过镜像检查、完整编译、单节点及双节点 smoke；
- 编译和运行审计能够证明二进制链接的是 SIF 内 OpenMPI；
- 文档准确写明 Slurm/PMI2、x86_64 和基础 TCP 运行边界；
- 任一不支持条件会前置失败，不会自动生成第二个站点镜像；
- pilogin 验证结果单独记录，但不进入镜像或 runner 命名。

## 11. 实施风险与处理

1. **镜像内 Slurm PMI2 client 与不同 Slurm 版本的兼容性**
   - 用跨 Slurm 主版本测试门禁验证；
   - 失败时调整唯一 PMI2 client pin，不能增加站点镜像。
2. **宿主模块环境污染容器**
   - 选择性保留 Slurm/PMI 变量；
   - 清除宿主 MPI/编译器路径，并用 `ldd`/`readelf` 强制审计。
3. **高速网络性能不一致**
   - 第一版以 TCP/shared memory 保证正确性；
   - UCX/OFI/RDMA 作为 runner 的显式、经测试优化，不作为默认保证。
4. **SIF 构建方式造成不可追溯差异**
   - 所有重建都从同一 OCI digest；
   - 发布 SIF SHA-256、构建日志和 SBOM。
5. **用户把“单一 SIF”误解为所有架构/调度器无条件运行**
   - README 明示受支持的架构和 launcher contract；
   - 其他 scheduler 通过 adapter 扩展，绝不通过复制镜像扩展。

## 12. 官方依据

- [OpenMPI 4.1 stable releases](https://www.open-mpi.org/software/ompi/v4.1/)
- [Slurm MPI Users Guide：OpenMPI 通过 `--with-pmi` 构建后可由 `srun --mpi=pmi2` 直接启动](https://slurm.schedmd.com/mpi_guide.html)
- [SchedMD Slurm 官方源码仓库](https://github.com/SchedMD/slurm)
