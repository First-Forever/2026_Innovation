# 2026 浙江大学“蒲公英”创新大赛——跨癌智筛

**Last updated：2026.03.15**

负责人、代码作者：李牧轩

项目成员：方颖博、向颖喆、陆梓萌、叶依朵、周昕悦

中文 | [English](#english)

## 项目简介
本仓库为大学生创新项目的代码与示例数据仓库，主题为基于转录组表达数据的恶性状态识别与机器学习预测。

当前仓库主要包含三部分内容：
- **R 预处理脚本**：对原始 RNA-seq count 数据进行 VST（variance stabilizing transformation）标准化；
- **Python 机器学习流程**：完成数据读取、特征筛选、GBM 模型训练、模型保存与测试；
- **示例测试数据**：提供 `TCGA-LAML` 数据集和 `meta.csv`，用于演示模型预测流程。

本仓库更适合作为项目原型、方法展示和流程复现示例，不直接作为临床诊断工具使用。

---

## 仓库结构

```text
2026_Innovation/
├── Data/
│   ├── TCGA-LAML-VST.txt.gz            # VST 标准化后的 LAML 测试数据
│   ├── TCGA-LAML.star_counts.tsv.gz    # 原始 LAML count 数据
│   └── meta.csv                        # 训练样本元数据
├── GBM_vst.pkl                         # 已训练好的 GBM 模型
├── ML_pre.ipynb                        # Python notebook：训练与测试流程
├── ML_pre.html                         # notebook 导出版本
└── VST.R                              # R 脚本：count 数据的 VST 标准化
```

---

## 方法流程

### 1. 数据预处理（R）
`VST.R` 使用 `DESeq2` 对 RNA-seq count 表达矩阵进行方差稳定化变换（VST），得到后续机器学习建模可直接使用的表达矩阵。

主要步骤包括：
1. 读取原始 count 矩阵；
2. 去除基因 ID 版本号；
3. 读取并对齐样本注释信息；
4. 基于 `DESeq2` 进行 VST 标准化；
5. 输出处理后的表达矩阵。

### 2. 模型训练（Python）
`ML_pre.ipynb` 中构建了基于 `scikit-learn` 的机器学习流程，核心模型为 **Gradient Boosting Classifier (GBM)**。

训练流程包括：
1. 读取 VST 标准化后的表达矩阵；
2. 读取 `meta.csv` 中的样本标签信息；
3. 以 `cancer_status` 作为预测目标；
4. 通过 `Pipeline` 串联特征过滤与模型训练；
5. 保存训练好的模型。

Notebook 中使用的主要策略：
- `VarianceThreshold(threshold=0.2)`：去除低方差特征；
- `SelectKBest(f_classif, k=800)`：筛选 800 个差异显著特征；
- `GradientBoostingClassifier`：进行分类建模。

### 3. 模型测试（Python）
仓库中提供了 `TCGA-LAML` 示例测试数据。`ML_pre.ipynb` 中加载 `GBM_vst.pkl` 后，对测试数据进行特征对齐与预测。

当前示例测试的目标是验证模型能够识别 AML 样本为 `Malignant`。

---

## 运行环境

### R 依赖
请确保已安装以下 R 包：
- `tidyverse`
- `data.table`
- `DESeq2`

### Python 依赖
建议 Python >= 3.9，并安装以下库：
- `pandas`
- `numpy`
- `matplotlib`
- `seaborn`
- `joblib`
- `polars`
- `scikit-learn`
- `jupyter`

可参考：

```bash
pip install pandas numpy matplotlib seaborn joblib polars scikit-learn jupyter
```

---

## 使用说明

### A. 运行 VST 预处理
若你有原始 count 数据，可先运行：

```r
source("VST.R")
```

该脚本会将原始 count 数据转换为 VST 标准化表达矩阵。

> 注意：`VST.R` 中目前包含本地工作目录设置和部分输入文件名，需要根据你的实际环境修改路径。

### B. 运行机器学习流程
打开 notebook：

```bash
jupyter notebook ML_pre.ipynb
```

然后依次执行：
1. 加载表达矩阵与样本信息；
2. 训练 GBM 模型；
3. 保存模型；
4. 加载 `GBM_vst.pkl`；
5. 对 `TCGA-LAML` 示例数据进行预测。

---

## 数据说明

### `meta.csv`
该文件包含训练样本的元数据，当前主要字段包括：
- `sample_id`：样本编号
- `subject_id`：受试者编号
- `sex`：性别
- `age_at_diagnosis`：诊断年龄
- `disease`：疾病类型
- `tissue`：组织来源
- `study`：研究来源
- `cancer_status`：标签（如 `Malignant`）

当前仓库中的 `meta.csv` 共包含 **12,352** 条样本记录，可用于建模标签与分组信息管理。

---

## 注意事项
- 本仓库主要用于科研训练、方法展示和项目申报展示；
- 仓库中部分脚本仍保留本地路径或实验性写法，复现前请根据实际情况调整；
- Notebook 中引用的训练输入文件名与仓库现有示例文件可能不完全一致，运行前建议先检查数据路径与文件名；
- 模型结果仅作为机器学习实验示例，不能替代真实临床判断。

---

## 后续可改进方向
- 增加 `requirements.txt` 或 `environment.yml`；
- 补充训练数据来源与预处理说明；
- 增加更多外部验证结果与性能评估图；
- 提供更清晰的输入/输出数据格式说明。

---

## English

## Overview
This repository contains code and example data for a student innovation project on malignant state recognition using transcriptomic expression profiles and machine learning.

The repository currently includes three main components:
- **R preprocessing script** for VST (variance stabilizing transformation) of raw RNA-seq count data;
- **Python machine learning workflow** for feature selection, GBM model training, model saving, and testing;
- **Example test data** including `TCGA-LAML` and `meta.csv` for demonstrating the prediction pipeline.

This repository is intended as a project prototype and workflow demonstration. It is **not** intended for direct clinical use.

---

## Repository Structure

```text
2026_Innovation/
├── Data/
│   ├── TCGA-LAML-VST.txt.gz
│   ├── TCGA-LAML.star_counts.tsv.gz
│   └── meta.csv
├── GBM_vst.pkl
├── ML_pre.ipynb
├── ML_pre.html
└── VST.R
```

---

## Workflow

### 1. Data preprocessing (R)
`VST.R` uses `DESeq2` to perform variance stabilizing transformation (VST) on RNA-seq count matrices.

Main steps:
1. Read raw count matrix;
2. Remove version suffixes from gene IDs;
3. Read and align sample metadata;
4. Perform VST normalization with `DESeq2`;
5. Export transformed expression matrix.

### 2. Model training (Python)
`ML_pre.ipynb` implements a machine learning workflow based on **Gradient Boosting Classifier (GBM)** using `scikit-learn`.

Training steps:
1. Load VST-normalized expression matrix;
2. Load labels from `meta.csv`;
3. Use `cancer_status` as the prediction target;
4. Train the model with a `Pipeline`;
5. Save the trained model.

Key methods used in the notebook:
- `VarianceThreshold(threshold=0.2)` for low-variance feature filtering;
- `SelectKBest(f_classif, k=800)` for selecting top informative genes;
- `GradientBoostingClassifier` for classification.

### 3. Model testing (Python)
The repository provides `TCGA-LAML` as an example external dataset. The notebook loads `GBM_vst.pkl`, aligns the features, and performs prediction on the test set.

In the current example, the AML samples are predicted as `Malignant`.

---

## Requirements

### R packages
- `tidyverse`
- `data.table`
- `DESeq2`

### Python packages
Recommended Python >= 3.9.
- `pandas`
- `numpy`
- `matplotlib`
- `seaborn`
- `joblib`
- `polars`
- `scikit-learn`
- `jupyter`

Install Python dependencies with:

```bash
pip install pandas numpy matplotlib seaborn joblib polars scikit-learn jupyter
```

---

## How to Use

### A. Run VST preprocessing
If you have raw count data, run:

```r
source("VST.R")
```

This script converts raw count data into a VST-transformed expression matrix.

> Note: `VST.R` currently contains a local working directory and file paths that should be adjusted before reuse.

### B. Run the machine learning workflow
Launch Jupyter Notebook:

```bash
jupyter notebook ML_pre.ipynb
```

Then execute the notebook step by step to:
1. Load expression matrix and metadata;
2. Train the GBM model;
3. Save the model;
4. Load `GBM_vst.pkl`;
5. Predict on the `TCGA-LAML` example dataset.

---

## Data Description

### `meta.csv`
This file contains metadata for the training samples. Major columns include:
- `sample_id`
- `subject_id`
- `sex`
- `age_at_diagnosis`
- `disease`
- `tissue`
- `study`
- `cancer_status`

The current `meta.csv` file contains **12,352** sample records.

---

## Notes
- This repository is mainly for research training, method demonstration, and project presentation;
- Some scripts still contain local paths or experimental settings and may need adjustment before reuse;
- The file names referenced in the notebook may not fully match the example files currently included in the repository;
- The model outputs are for machine learning demonstration only and should not be interpreted as clinical diagnosis.

---

## Possible Future Improvements
- Add `requirements.txt` or `environment.yml`;
- Document the source and preprocessing of training data more clearly;
- Add more external validation and performance metrics;
- Provide a clearer specification of input/output formats.
