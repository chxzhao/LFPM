# LFPM：论文方法与复现代码

依据 `LFPM_revised_20260929` 整合稿整理。作者：Chenxi Zhao。
维护者邮箱：chxzhao@mail.nankai.edu.cn。

- `lfpm_fit()`：对收益面板拟合带截距的因子回归。
- `l_test()`：指定阶数的 L 检验。
- `adal_test()`：完整候选网格的 AdaL；显式选择模拟或实证校准。
- `lfpm_test()`：模拟章节口径下的 AdaL、PY、MAX、COM、max-k。
- `reproduce_simulation()`：短流程检查或完整论文模拟，可恢复计算。
- `reproduction_reference()`：读取论文模拟的汇总参考表；不声称重新模拟。
- `lfpm_empirical_test()`：通过保留原随机数与数值口径的 Python 后端，分析外部收益与因子矩阵。
- `reproduce_empirical()`：读取包外数据并执行论文实证流程。
- `report_empirical()`：从新计算的实证结果生成论文汇总表和 alpha 分布图。

在本目录运行 `R CMD build .`，再运行 `R CMD INSTALL LFPM_0.1.0.tar.gz`。
核心方法和模拟只依赖 R 标准包。实证的可选 Python 依赖、编译要求、数据来源、
输入格式和示例详见 `inst/reproduction/`。原始及清洗后的市场数据都不随包分发。
中国 CSMAR 数据需相应访问权限，数据源更新后也不保证与论文历史文件逐位一致。

模拟与实证的 AdaL 分量概率和端点截断存在明确差异，本包按整合稿保留二者。
不要用后来全经验尾的模拟版本替换该稿数字。其他比较方法大多是本项目的实现或
改编，不会写成调用不存在的作者官方 R 包；QUIC 等第三方来源与许可证单独保留。

正式模拟预算：每点评估 1000 次，独立零假设校准 1999 次，每次 bootstrap 999 次。
`smoke` 只用于安装与流程检查，不等于完成论文复现。完整说明见英文 README。

已知历史缺口：最终中国 max-k 惩罚参数敏感性的归档结果缺少原始随机种子和
预算记录。复现函数会以明确记录的种子真实重算，不把历史小数填入新结果。
