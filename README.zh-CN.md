# LFPM：论文方法、模拟与实证复现

对应论文 *Testing factor pricing restrictions in high dimensions without knowing the sparsity of alphas*。
版本 **0.1.2**。作者：Chenxi Zhao。维护者邮箱：chxzhao@mail.nankai.edu.cn。

本包提供可直接调用的 R 检验函数、论文模拟代码和实证复现代码。
**实际数据不放在包内**；数据来源链接、所需文件及字段约定见
[DATA_SOURCES.md](inst/reproduction/DATA_SOURCES.md)。数据和运行结果应存放在包及公开仓库之外。

## 安装与依赖

在本目录运行：

```sh
R CMD build .
R CMD INSTALL LFPM_0.1.2.tar.gz
```

需要 R 4.2 或以上版本。通用 R 方法和模拟仅使用 R 标准包。
实证复现另需 Python 3.9+，以及 NumPy、pandas、SciPy、statsmodels、threadpoolctl；
绘图需 matplotlib。实证 max-k 的 QUIC 求解器在运行时编译，需要 C++ 编译器，
macOS 使用 Accelerate，Linux 使用 LP64 BLAS/LAPACK。该求解器目前不支持 Windows，
可设置 `include_maxk=FALSE` 明确省略该比较方法。安装 R 包本身不编译可选求解器。
环境配置、版本记录和具体命令见[实证指南](inst/reproduction/EMPIRICAL.md)。

## 方法函数

- `simulate_lfpm()`：生成指定模型、创新分布和稀疏度的合成面板。
- `lfpm_fit()`：拟合带截距的因子回归。
- `l_test()`：指定阶数的 L 检验。
- `adal_test()`：完整候选网格的 AdaL，可选择分量校准方式。
- `lfpm_test()`：模拟章节口径下的 AdaL、PY、MAX、COM、max-k。
- `lfpm_empirical_test()`：在用户提供的收益和因子矩阵上调用论文实证实现。

```r
library(LFPM)
dat <- simulate_lfpm(model = 1, N = 40, T_size = 60,
                     distribution = "Normal", m = 0, seed = 42)
adal_test(dat$Y, dat$F, B = 99, calibration = "simulation", seed = 43)
adal_test(dat$Y, dat$F, B = 99, calibration = "empirical", seed = 43)
```

以上使用合成数据。`"simulation"` 保留论文模拟的混合分量校准：固定阶数 5
使用经验 bootstrap 尾概率，其余分量使用 bootstrap 均值和标准差下的正态尾。
`"empirical"` 表示所有分量使用经验 bootstrap 尾概率；它也是通用 R 算法选项。
两者的 Cauchy 端点截断不同。论文实证流程则保留 Python/NumPy 的原随机流；
R 与 NumPy 使用同一个整数种子不会产生同样的随机抽样。示例 `B=99` 仅用于快速演示。

## 模拟复现

- `simulation_design()`：查看或选择模拟场景。
- `reproduce_simulation()`：短流程检查或完整论文模拟，支持检查点恢复。
- `plot_simulation()`：绘制模拟结果。
- `reproduction_reference()`：读取缓存的论文模拟汇总表，不重新模拟。

```r
smoke <- reproduce_simulation("reproduction-output/smoke", preset = "smoke")
# 完整预算，计算量较大：
# paper <- reproduce_simulation("reproduction-output/paper", preset = "paper")
```

完整设计包含 96 个 size 场景、3200 个 power 场景，并保留全部 `m=1:200`。
正式预算为每点评估 1000 次、独立零假设校准 1999 次、每次 bootstrap 999 次。
`smoke` 只检查安装与流程，不等于完成正式模拟。
细节见[模拟指南](inst/reproduction/SIMULATION.md)。

## 实证复现

先从[数据来源页](inst/reproduction/DATA_SOURCES.md)获取外部数据，
按该页说明提供文件及字段。CSMAR 需相应访问权限；公开来源会更新，
当前下载文件不保证与论文历史版本完全一致。

```r
py <- "/path/to/your/python"
china <- reproduce_empirical(
  data_dir = "/external/market-data",
  output_dir = "/external/results/china",
  analysis = "china", python = py)
report <- report_empirical(china$output_dir)
```

`analysis="all"` 覆盖已维护的最终中国、美国 FF3/FF5 设置，以及诊断、
拟合零假设校准、断点重选和 max-k 惩罚参数敏感性计算，计算量较大。
主设置运行后会自动生成联合检验、分量概率、alpha 经济量级、BH 数量、
已有诊断和分布图。`report_empirical()` 可从已计算结果重建报告。
缺失市场或省略方法会明确标识，不以历史论文数值补齐。
各参数、预算、种子和输出见[实证指南](inst/reproduction/EMPIRICAL.md)。

## 方法归属与复现边界

比较方法是本项目的实现或改编，不会声称调用作者官方 R 包。
`cvCovEst::thresholdingEst` 仅是原项目可选数值核验参照，并非运行依赖。
实证 max-k 使用随包保留来源及许可证的 QUIC 1.2 源码；统计文献和软件归属见
[方法说明](inst/reproduction/METHODS.md)。

短流程和合成数据测试不能证明完整模拟及真实数据分析已经重新运行，
也不代表所有论文数值都已重新验证。最终中国 max-k 惩罚参数敏感性的一项
历史结果缺少原始种子和抽样预算记录；代码提供明确种子下的真实重算，
标记 `historical_rng_parity_verified=FALSE`，不承诺完全还原那个 Monte Carlo 小数。
这一缺口与已恢复种子规则的主检验流程分开说明。

包内不含原始数据、清洗面板、数据库导出或凭据；运行输出也可能含个股标识和
估计量，请存放在公开仓库之外。许可证：GPL-3，QUIC 保留其上游许可证说明。
