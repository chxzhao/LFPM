# LFPM：论文方法与模拟复现

版本 **0.1.1**，依据 `LFPM_revised_20260929` 整合稿整理。
作者：Chenxi Zhao。维护者邮箱：chxzhao@mail.nankai.edu.cn。

本包提供可直接调用的 R 检验函数和论文模拟复现。
实际案例部分只提供[数据来源链接](inst/reproduction/DATA_SOURCES.md)，
不包含实际案例分析代码、复现或报告函数，也不包含原始或清洗后的市场数据。

## 安装

在本目录运行：

```sh
R CMD build .
R CMD INSTALL LFPM_0.1.1.tar.gz
```

需要 R 4.2 或以上版本，仅使用 R 标准包，不需要编译扩展。

## 通用方法函数

- `simulate_lfpm()`：生成指定模型、创新分布和稀疏度的合成面板。
- `lfpm_fit()`：拟合带截距的因子回归。
- `l_test()`：指定阶数的 L 检验。
- `adal_test()`：完整候选网格的 AdaL，可明确选择分量校准方式。
- `lfpm_test()`：模拟章节口径下的 AdaL、PY、MAX、COM、max-k。

```r
library(LFPM)
dat <- simulate_lfpm(model = 1, N = 40, T_size = 60,
                     distribution = "Normal", m = 0, seed = 42)
adal_test(dat$Y, dat$F, B = 99, calibration = "simulation", seed = 43)
adal_test(dat$Y, dat$F, B = 99, calibration = "empirical", seed = 43)
```

以上使用合成数据。`"simulation"` 保留论文模拟的混合分量校准：固定阶数
5 使用经验 bootstrap 尾概率，其余分量使用 bootstrap 均值和标准差下的正态尾。
`"empirical"` 表示所有分量使用经验 bootstrap 尾概率，是通用算法选项，
不代表本包提供实际案例分析流程。两者的 Cauchy 端点截断也不同，详见英文 README
和函数帮助。示例 `B=99` 只用于快速演示。

## 模拟复现

- `simulation_design()`：查看或选择模拟场景。
- `reproduce_simulation()`：运行短流程检查或完整论文模拟，支持检查点恢复。
- `plot_simulation()`：绘制模拟结果。
- `reproduction_reference()`：读取缓存的论文模拟汇总参考表；不会重新模拟。

```r
smoke <- reproduce_simulation("reproduction-output/smoke", preset = "smoke")
# 完整预算，计算量较大：
# paper <- reproduce_simulation("reproduction-output/paper", preset = "paper")
```

完整设计包含 96 个 size 场景、3200 个 power 场景，并保留全部 `m=1:200`。
正式预算为每点评估 1000 次、独立零假设校准 1999 次、每次 bootstrap 999 次。
`smoke` 只用于安装与流程检查，不等于完成正式模拟。
具体场景选择、恢复规则和输出说明见[模拟指南](inst/reproduction/SIMULATION.md)。

其他比较方法是本项目的实现或改编，不会声称调用作者官方 R 包。
原项目曾用 `cvCovEst::thresholdingEst` 做可选数值核验，它不是运行依赖。
方法出处见[方法归属说明](inst/reproduction/METHODS.md)。

数据来源由各提供方维护，访问和使用需遵守相应条件；案例数据和代码不随包分发。
许可证：GPL-3。
