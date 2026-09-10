# Fourier_Quad C++ 流水线指南

本文说明如何编译、配置和运行 `Fourier_Quad_Cpp`。

> English: [CPP_GUIDE.md](CPP_GUIDE.md)

## 程序结构

[`cpp_Standard`](cpp_Standard/) 与 [`cpp_Lite`](cpp_Lite/) 均生成
`Fourier_Quad_Pipe`，并包含：

- `main.cpp`：MPI 生命周期与六个顶层阶段的调度；
- `pipeline.example.ini`：完整运行期配置模板；
- `config/`：编译默认值与固定数值参数；
- `include/`、`src/`：阶段模块、运行配置与星表布局；
- `include/general/`、`src/general/`：曝光表、路径、MPI、调度与通用数值工具；
- `Makefile`：C++17/MPI 构建入口。

Standard 保留平场、掩膜、简化测天、外部/混合/PCA PSF 等可选分支。Lite 只保留
Gaia 测天、逐 CCD DQ 掩膜、外部源星表、帧内恒星、去混叠、局域多项式 PSF，且
不含 PCA。Lite 是物理删除未使用分支后的版本，不是换了默认值的 Standard。

## 处理流程

六个顶层阶段按固定顺序执行：

| 阶段 | CLI | 作用 |
|---|---|---|
| `process_astrocat` | `--run-astrocat` | 将原始两列 Gaia 星表重分块、去重为一度瓦片。 |
| `process_extcat` | `--run-extcat` | 将原始 External source catalog 文件重分块为天空瓦片。 |
| `process_init` | `--run-init` | 搜索 `.fits.fz`，提取 Science/DQ CCD，生成曝光表。 |
| `process_main` | `--run-main` | 执行九阶段数值流水线。 |
| `process_rearr` | `--run-rearr` | 将 `*_all.cat` 空间分区。 |
| `process_fd` | `--run-fd` | 按场畸变分箱恢复平均剪切。 |

`process_astrocat`、`process_extcat` 每次运行依次各执行一次，其余阶段按数据集顺序执行；
首次集体失败会停止整个任务。

`process_main` 用素数乘积选择阶段：

| 阶段 | 素数 | 内容 |
|---:|---:|---|
| 1 | 2 | 背景/噪声预处理与 Gaia 匹配 |
| 2 | 3 | 测天解 |
| 3 | 5 | 源检测、去混叠、恒星候选体 |
| 4 | 7 | 恒星候选体功率谱 |
| 5 | 11 | PSF 选择与建模 |
| 6 | 13 | 星系功率谱 |
| 7 | 17 | Fourier_Quad 估计量与形态测量 |
| 8 | 19 | 曝光级统计 |
| 9 | 23 | 星表合并与标定 |

`process_stage` 能被某素数整除时执行对应阶段。完整值 `223092870` 启用全部阶段；
阶段 9 必须与阶段 8 同时启用。

## 编译

需要支持 C++17 的 MPI C++ 编译器，以及 CFITSIO、FFTW3、Eigen3、LAPACK、BLAS。

```bash
cd cpp_Lite                    # 或 cpp_Standard
make -j4
./Fourier_Quad_Pipe --help
```

科学库集中安装在同一前缀时：

```bash
make CXX=/path/to/mpicxx STACK_PREFIX=/opt/science-stack -j4
```

Eigen 单独安装时再传入 `EIGEN_INCLUDE=/opt/eigen/include/eigen3`。当前 Standard 与
Lite Makefile 只公开 `all` 和 `clean`；`make` 等同于 `make all`。修改编译期配置或
工具链后，先执行 `make clean` 再重新编译。

## 配置

支持3种配置方式：
1. 配置文件 `pipeline.ini`；

复制所选版本的模板，修改为所需配置：

```bash
cp pipeline.example.ini pipeline.ini
```

2. 命令行参数；
3. `config`下头文件中设置编译器常量。

优先级为：

```text
编译默认值 < INI < CLI
```

INI 分为五段：`[process]` 控制阶段和输出路径，`[astrocat]` 控制原始两列 Gaia 输入、
独立输出目录、header 与已有输出策略，`[extcat]` 控制原始星表解析与布局，`[init]`
控制归档根目录和数据集，`[lensing]` 控制可在运行间改变的科学分支与相机参数。
Lite 会拒绝 Standard 专属键；未知段、未知键、无效值及不一致的阶段/schema 组合均会在
执行前报错。

CLI 同时接受 `--name value` 与 `--name=value`。布尔值接受 `true/false`、`1/0`、
`yes/no`、`on/off`。首次显式 `--dataset`、`--contains`、`--extcat-contains` 会替换
INI 列表，后续同名选项追加。一个裸位置参数仍可作为 `--expo-list` 的兼容写法。

完整 CLI 以 `./Fourier_Quad_Pipe --help` 为准；所有 INI 键见
[CPP_PIPELINE_PARAMETERS.md](CPP_PIPELINE_PARAMETERS.md)。

运行期设置优先写入所选版本的 INI，CLI 只用于单次覆盖。头文件专属设置应修改
`cpp_Standard/config/` 或 `cpp_Lite/config/` 下对应的配置头文件；修改源参数而非
派生常量，然后执行 `make clean` 并重新编译。参数参考按配置领域分别给出完整的
Standard/Lite 对照表，并明确标出哪些值可在不重编译的情况下覆盖。

### 集中路径配置

所选版本的 `config/pathconfig.hpp` 是所有输入/输出路径、流程输出与曝光表名称、重排
文件名，以及初始化/处理产物相对目录的唯一实际定义点。原有命名空间与 RuntimeConfig
键保持不变。

- 必须保留 `AstroCatConfig::ASTROCAT_OUTPUT_DIRECTORY` 从
  `LensingConfig::ASTROMETRY_CAT` 初始化的关系，以及这两个相互耦合但独立的符号。
  RuntimeConfig 构造后，
  `[astrocat].output_directory` 与 `[lensing].astrometry_cat` 是相互独立的生产者/消费者
  字段；若它们应指向同一批新瓦片，必须同时显式配置。
- 必须保留 `ExtCatConfig::EXTCAT_OUTPUT_DIRECTORY` 对
  `LensingConfig::SOURCE_CAT` 的 `const std::string&` 引用关系。
  `[lensing].source_cat` 与 `[extcat].output_directory` 更新同一个有效 RuntimeConfig 目录，
  `--extcat-output` 具有最终优先级。
- `FLAT_PATH` 与 `PSF_PATH` 只存在于 Standard；Lite 已物理删除对应可选分支。
- `ProcessRearrConfig` 的固定文件名和两组 `OutputLayout` 目录数组没有
  RuntimeConfig/INI/CLI 覆盖。直接修改它们或 `pathconfig.hpp` 中任何编译默认值后，
  必须执行 `make clean && make`。

INI 与 CLI 只覆盖 RuntimeConfig 副本，不会修改头文件或其编译期引用关系。

### 常用及随图像数据源变化的参数

下表是运行前应主动检查的参数入口。运行时设置优先写入 INI；有对应 CLI 的字段可在
单次调用中覆盖。只有标为“编译时”的固定参数才需要修改所选版本的头文件并重新执行
`make clean && make`。派生尺寸和运行时解析出的星表列号不要单独修改。

| 类别 | 参数（当前默认） | 修改方式 | 何时修改与约束 |
|---|---|---|---|
| 顶层阶段 | `[process].run_process_astrocat/run_process_extcat/run_process_init/run_process_main/run_process_rearr/run_process_fd` | INI；运行时 `--run-astrocat`、`--run-extcat`、`--run-init`、`--run-main`、`--run-rearr`、`--run-fd` | 选择本次执行的阶段。Standard 默认 `false/false/true/true/true/true`，Lite 默认 `false/false/true/true/false/false`。 |
| Science/DQ 归档与数据集 | `[init].science_root`、`dq_root`、`output_root`、`datasets`、`contains` | INI；运行时 `--science-root`、`--dq-root`、`--output-root`、`--dataset`、`--contains` | 更换观测归档、文件名前缀、筛选 token 或输出根目录时修改。Lite 必须提供逐 CCD DQ masks。 |
| 曝光表与阶段输出 | `[process].expo_list`、`rearr_output_directory`、`rearr_output_base_directory`、`rearranged_expo_list_filename`、`rearranged_expo_list_directory`、`fd_expo_list`、`fd_output_directory`、`fd_output_base_directory` | INI；对应 CLI 为 `--expo-list`、`--rearr-output-dir`、`--rearr-output-base`、`--rearr-list-name`、`--rearr-list-dir`、`--fd-expo-list`、`--fd-output-dir`、`--fd-output-base` | 下游单独运行，或改变重排/FD 输出目录和曝光表位置时修改。 |
| 固定生成布局 | `SKIP_DIRECTORY_NAME`、`SUBCAT_PREFIX`、`SUBCAT_EXTENSION`、`SUMMARY_FILENAME`、`NON_CHIP_BASE_DIRECTORIES`、`CHIP_PRODUCT_DIRECTORIES` | `config/pathconfig.hpp`，编译时 | 仅在发布星表命名或相对输出目录约定变化时修改；重编译并重新生成受影响产物。 |
| Gaia 星表分块 | `[astrocat].input_directory`、`output_directory`、`add_header=true`、`existing_policy=fail` | INI；运行时 `--astrocat-input`、`--astrocat-output`、`--astrocat-add-header`、`--astrocat-existing` | 更换 Gaia 原始星表或重跑策略时修改。输出只属于 `process_astrocat`，不会传播到 `[lensing].astrometry_cat`。 |
| Gaia 星表布局 | `[lensing].astrometry_cat_type=1`、`astrometry_cat`；`ASTROMETRY_TILE_PREFIX="astra_"` | 布局/路径由 INI 运行时设置；前缀在 `config/pathconfig.hpp` 编译时设置 | `1` 读取旧式大 `gaia_*.cat` 瓦片；`2` 累积读取 `process_astrocat` 生成的一度 `<prefix>RA_*.dat` 瓦片。切换布局时同时指向消费目录；修改前缀后需重编译。 |
| 外部星表发现与发布 | `[extcat].input_directory`、`output_directory`；`SOURCE_CAT_TILE_PREFIX="extern_"` | 路径由 INI/CLI 运行时设置；前缀在 `config/pathconfig.hpp` 编译时设置 | 更换 External source catalog 的原始目录或规范化瓦片目录时修改。输出目录不能等于或位于输入目录内；它同时是有效的 `SOURCE_CAT`，生产者与消费者共用该前缀。 |
| 外部星表 schema | `[extcat].total_columns`、`use_explicit_columns`、`input_columns`、`use_explicit_coordinate_columns`、`ra_column`、`dec_column`、`zp_column` | INI；投影和 RA/Dec/ZP 列可用 `--extcat-columns`、`--extcat-ra-column`、`--extcat-dec-column`、`--extcat-zp-column` 覆盖 | 更换 survey 或列顺序时修改。显式投影必须保留 RA、Dec、photo-z 以及启用阶段需要的字段；完整行宽和下游偏移由 `CatalogLayout` 自动解析。 |
| FD magnitude 映射 | `[extcat].mag_g_column`、`mag_r_column`、`mag_i_column`、`mag_z_column`、`mag_y_column` | INI，运行时；`0` 表示该波段不存在 | 更换 survey/band schema 时修改。显式投影应保留所用 magnitude；`process_fd` 至少需要一个波段，并按 `i -> z -> r -> g -> y` 选择首个可用列，无需手工修改 FD 偏移。 |
| 标定路径与 Standard 分支 | `[lensing].flat_path`、`psf_path`、`astrometry_trivial=0`、`include_flat=0`、`include_mask=2`、`ext_cat=1`、`ext_psf=0`、`psf_type=1`、`psf_ms=0` | Standard 的 INI 运行时设置；Lite 拒绝这些已删除分支键 | 更换平场/外部 PSF 数据或选择替代科学分支时修改。Lite 固定为 Gaia、无平场、逐 CCD DQ、外部源星表、帧内 PSF、局域多项式且无 PCA。 |
| 图像与探测器几何 | `[lensing].ccd_split=2`、`pixel_size=0.2628`、`nmax_chip=62`、`chipnx=2046`、`chipny=4094`；编译时 `ns=64`、`chip_margin=8` | 前五项为 INI 运行时；后两项在 `config/LensingConfig.hpp` 中修改并重编译 | 更换相机、放大器布局、像元尺度、PSF map 几何或 stamp 尺寸时成组核对。Stage 1 动态读取 Science FITS 轴长；成功的 Standard Hybrid PSF 图与 FD 边界使用有效运行时 `chipnx/chipny`。本仓库没有固定 `npx/npy`。 |
| 数值阶段 | `[lensing].process_stage=223092870` | INI，运行时 | 用素因数选择九个主流程阶段；阶段 9（23）必须与阶段 8（19）同时启用。 |
| 源检测与像素阈值 | `saturation_thresh=25000` | `config/LensingConfig.hpp`，编译时 | 更换图像源、增益或饱和定义后，以代表性数据重新标定并重编译。 |
| FD 探测器规则 | `bad_ccds={}`、`chip_mask_edge=50` | `config/FDConfig.hpp`，编译时 | 接受闭区间按有效运行时几何自动生成为 `[edge, chipnx-edge]` 与 `[edge, chipny-edge]`。更换相机、坏 CCD 清单或边缘策略时修改；`n_bad_ccds` 为派生长度。 |

每个独立参数的 Standard/Lite 默认值、合法值、INI/CLI 覆盖和重编译要求见
[CPP_PIPELINE_PARAMETERS.md](CPP_PIPELINE_PARAMETERS.md)。修改高耦合参数时，应保留基准
配置，并先用最小代表性数据验证。

## 常用运行方式

按 INI 串联运行：

```bash
mpirun -np 8 ./Fourier_Quad_Pipe --config pipeline.ini
```

仅对已有曝光表运行主流程：

```bash
mpirun -np 8 ./Fourier_Quad_Pipe --config pipeline.ini \
  --run-init false --run-main true --run-rearr false --run-fd false \
  --expo-list /data/work/expo_gband.list
```

仅重分块原始 Gaia 星表：

```bash
mpirun -np 8 ./Fourier_Quad_Pipe \
  --run-astrocat true --run-extcat false --run-init false --run-main false \
  --run-rearr false --run-fd false \
  --astrocat-input /data/raw_gaia --astrocat-output /data/gaia/tiles \
  --astrocat-add-header true --astrocat-existing fail
```

`--astrocat-output` 只控制 `process_astrocat`，不会与
`[lensing].astrometry_cat` 比较、同步，也不要求两者一致。

仅生成外部星表瓦片：

```bash
mpirun -np 8 ./Fourier_Quad_Pipe \
  --run-extcat true --run-init false --run-main false \
  --run-rearr false --run-fd false \
  --extcat-input /data/raw_catalogs --extcat-output /data/catalogs/tiles
```

至少启用一个顶层阶段。初始化成功后，后续阶段自动使用生成的绝对路径
`expo_<target>.list`。在 Slurm 上应改用站点支持的启动器；本仓库容器 runner 固定使用
`srun --mpi=pmi2`。

## 输入与输出

Science images、Gaia catalog、External source catalog 与 DQ masks 的统一定义和最低
要求见根目录 README 的[输入数据要求](README_CN.md#输入数据要求)。DQ masks 是取决于
配置的输入类别：Lite 始终读取逐 CCD DQ masks；Standard 只有在所选
`[lensing].include_mask` 模式和实际执行路径都不访问 DQ 数据时才可省略。

曝光表每个非空记录包含一个 CCD 列表路径，可带兼容用的 CCD 数量；当前 C++ 读取器
支持带引号的路径。

初始化器读取原始压缩归档但不移动它们，在 `output_root/<target>/` 下创建 `science/`、
`dqmask/`、`stamps/`、`expolists/`、`result/`，并发布 `expo_<target>.list`、
`fits_<target>.list` 与 manifest。逐曝光 CCD 列表保存在 `<dataset>/expolists/`，
数值阶段产物仍保存在 `<dataset>/stamps/`。主要下游结果为：

```text
<dataset>/result/<exposure>_all.cat
<dataset>/<rearr-output-dir>/subcat_*.cat
<dataset>/<fd-output-dir>/FD_test_comb.dat
```

`_all.cat`是按曝光为单位的剪切目录，默认包含外部星表字段、原始 1-based
`EXPO_NUM`、1 个 CCD 编号和 25 个流水线字段，共 45 列。schema 升级后需重新生成
Stage 9、rearr 与 FD 产物；旧 44 列数据不能与新版混用。
`subcat_*.cat`是按RA/DEC 重新分块的目录，单个源的所有测量记录连续排列，便于快速去重。
`FD_test_comb.dat`是程序process FD生成的场畸变测试表格文件，用于矫正剪切测量。

## 常见错误

- 不要在未启用阶段 8 时启用阶段 9。
- 外部星表输出目录不能等于或位于输入目录之下。
- 显式投影必须包含原始 RA、Dec、photo-z 列。
- Lite 不能使用 Standard 专属 `[lensing]` 键。
- 容器内程序参数必须使用容器路径，而不是宿主路径。
- 不要让多个作业同时清理或编译同一源码副本。

本地 Docker 见 [cpp_docker/README-CN.md](cpp_docker/README-CN.md)，Slurm/Apptainer
见 [cpp_docker/runner/README-CN.md](cpp_docker/runner/README-CN.md)。镜像只提供工具链，
不包含流水线源码、配置、星表或观测数据。
