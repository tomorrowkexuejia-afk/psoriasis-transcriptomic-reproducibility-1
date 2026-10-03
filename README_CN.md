# 银屑病转录组可重复性研究公开仓库说明（v1.0）

这是用于 Zenodo 等永久仓库上传的统一可重复性分析包。核心分析与 2026-09-17 新增的 paired audit-impact analysis 已整合在同一个目录中。

## 最重要的目录

- `core_analysis/code/`：分析与核验代码
- `core_analysis/inputs/`：处理后的输入矩阵、样本映射和冻结输入
- `core_analysis/evidence/`：GEO 元数据、样本重叠证据、来源判定材料
- `core_analysis/results/`：核心分析数值结果
- `core_analysis/tables/`：claim/cohort 级结果表
- `core_analysis/figures/`：主图及 Figure S1（PDF/PNG/SVG）
- `core_analysis/audit_impact/`：来源审计、样本去重、方向冻结和推断层级的配对比较
- `historical_notes/`：历史稿件和旧的稿件生成辅助文件，仅用于版本追溯，不代表当前投稿正文

## 当前主分析集

最终 source-adjudicated development synthesis 为：13 篇来源论文、33 条 paper-gene claims、32 个基因、176 个 cohort contributions。

GSE54456 与 GSE13355 共用的 42 个样本已从所有最终 eligible GSE54456 contributions 中删除，最终 GSE54456 每次使用 72 例 psoriasis + 60 例 healthy。

## 一键重跑

```bash
cd core_analysis
pip install -r requirements.txt
python code/run_pipeline.py
python code/audit_impact_analysis.py
```

也可以在根目录执行 `python run_all.py`。

正式上传 Zenodo 后，请把获得的 DOI 填入投稿稿件中的 Data availability / Code availability。
