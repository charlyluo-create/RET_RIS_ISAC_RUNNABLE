RET-RIS ISAC 完整可运行代码包
================================

一、主入口
在 MATLAB 切换到本目录后运行：

    run_all_paper_experiments

如中途停止，并且前面实验的 results/*.mat 已经存在，可从指定实验继续：

    run_all_paper_experiments(6)

二、运行环境
- MATLAB R2025b 或兼容版本
- CVX 已安装并执行过 cvx_setup
- Parallel Computing Toolbox 可选；没有并行工具箱时实验脚本会按自身逻辑退回串行或提示

三、默认流程
主入口依次完成：
1. 最小依赖检查
2. 通信功率求解器 duality 与 CVX 交叉验证
3. Experiment 1--12
4. MA-RCG-MF V2 运行后审计
5. 导出论文统计表
6. 生成 6 张组合实验图

四、图片格式
本版本的 generate_paper_figures.m 按用户提供的参考格式生成：
- Figure 1：2×2，共 4 个子图
- Figure 2：2×2，共 4 个子图
- Figure 3：2×2，共 4 个子图
- Figure 4：2×2，共 4 个子图
- Figure 5：1×2，共 2 个子图
- Figure 6：2×2，共 4 个子图
合计 6 张组合图、22 个子图。

统一特点：
- 白底；
- Times New Roman；
- 分类横轴采用倾斜标签以减少重叠；
- 多方案曲线保持一致的颜色、线型与 marker；
- 误差条和置信区间保留；
- PNG 600 dpi；
- PDF 矢量输出；
- 同时保存可编辑 MATLAB FIG。

生成目录格式：

    Main_Figures_paper_yyyymmdd_HHMMSS/

五、核心文件
- run_all_paper_experiments.m：唯一推荐总入口
- unified_ret_ris_isac_config.m：统一系统/实验参数
- get_figure_policy.m：单实验绘图策略
- generate_paper_figures.m：6 张组合图生成器
- export_all_paper_statistics.m：统计表导出
- exp1_* 到 exp12_*：12 组实验
- function/：所有实际依赖的公共函数

六、说明
本包已统一移除期刊缩写相关的文件名、函数名、输出目录名、注释和命令行提示，避免后续改投其他期刊时再修改代码。
实验模型、MA-RCG-MF V2 算法、strict-paired mobility 逻辑和数值参数保持不变；本次主要统一工程命名与最终图片生成格式。
