# DimerPlus - 引物二聚体分析工具

一个基于 Streamlit 的引物二聚体分析和优化工具，支持批量序列分析、智能合孔方案生成等功能。

## 功能特点

- 📊 批量引物二聚体分析
- 🧬 基于 Primer3 的热力学计算
- 🔍 智能筛选和排序
- 🎯 OR-Tools 全局优化算法
- 📋 多种合孔方案生成
- 📥 结果导出（Excel/CSV）

## 安装依赖

### 本地运行

```bash
pip install -r requirements.txt
```

### 主要依赖

- `streamlit` - Web 应用框架
- `pandas` - 数据处理
- `biopython` - 生物信息学工具
- `primer3-py` - 引物设计和热力学计算
- `ortools` - Google 优化工具（全局最优算法）
- `openpyxl` - Excel 文件处理

## 运行应用

```bash
streamlit run dimerplus.py
```

## Streamlit Cloud 部署

1. 将代码推送到 GitHub 仓库
2. 确保仓库包含 `requirements.txt` 文件
3. 在 [Streamlit Cloud](https://share.streamlit.io/) 创建新应用
4. 选择您的 GitHub 仓库和主文件 `dimerplus.py`
5. 部署完成后即可使用

## 使用说明

1. **上传序列文件**：支持 Excel 格式，需包含序列名称和序列信息
2. **设置分析参数**：配置温度、盐浓度等条件
3. **运行二聚体分析**：自动计算所有序列对的二聚体形成情况
4. **筛选和排序**：根据 ΔG 值、Tm 值等条件筛选结果
5. **生成合孔方案**：使用智能算法生成最优的引物组合方案
6. **导出结果**：下载分析结果和合孔方案

## 故障排除

### OR-Tools 安装问题

如果遇到 "请安装Google OR-Tools" 错误：

**本地环境：**
```bash
pip install ortools
```

**Streamlit Cloud：**
- 确保 `requirements.txt` 包含 `ortools>=9.7.0`
- 重新部署应用

**其他云平台：**
- 检查部署配置是否包含所有依赖
- 确保 Python 版本兼容（推荐 3.8+）

### 常见问题

1. **文件上传失败**：检查文件格式是否为 Excel，且包含正确的列名
2. **计算超时**：减少序列数量或调整参数设置
3. **内存不足**：在云平台上可能需要升级资源配置

## 技术架构

- **前端**：Streamlit Web 界面
- **计算引擎**：Primer3-py（热力学计算）
- **优化算法**：Google OR-Tools CP-SAT 求解器
- **数据处理**：Pandas + BioPython

## 许可证

本项目仅供学术研究使用。
