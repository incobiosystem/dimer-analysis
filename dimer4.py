import streamlit as st
import pandas as pd
import primer3
from Bio import Seq
import os

# 设置页面标题和图标
st.set_page_config(page_title="Dimer Analysis", page_icon="🔬", layout="wide")  # 使用宽布局

# 顶部布局，三列布局，带有间隔
col1, col_spacer1, col2, col_spacer2, col3 = st.columns([1, 0.5, 3, 0.5, 5])  # 使用间隔列来调整间距

# 左侧导航栏
with col1:
    logo_path = "logo.png"  # 确保logo文件名与路径正确
    if os.path.exists(logo_path):
        st.image(logo_path, width=300)  # 调整图片宽度
    st.markdown("<h1 style='font-size: 56px;'>Dimer Analysis</h1>", unsafe_allow_html=True)  # 可调整标题字体大小
    st.markdown("<h1 style='font-size: 24px;'>Incobiosystem</h5>", unsafe_allow_html=True)  # 可调整标题字体大小


# 中间列：上传文件和选择序列
with col2:
    st.subheader("上传需要分析的Excel文件:")
    uploaded_file = st.file_uploader("选择文件", type=["xlsx"])
    # 输入 Tm 值
    tmin = st.number_input("本次分析将输出高于以下温度的dimer (默认10℃):", value=10.0)

    if uploaded_file is not None:
        df = pd.read_excel(uploaded_file)
        sequences = []

        # 显示序列列表并添加多选框
        st.subheader("选择序列进行分析:")
        for index, row in df.iterrows():
            if "P" in row['物料名称']:
                selected = st.checkbox(row['物料名称'])
                if selected:
                    sequences.append(row['物料名称'])

    

        
with col3:
# 开始分析按钮
    if st.button("开始分析"):  # 可以调整按钮大小和样式
        if len(sequences) < 2:
            st.warning("请至少选择两个序列进行分析！")
        else:
            dcP = {row['物料名称']: row['序列'] for index, row in df.iterrows() if row['物料名称'] in sequences}
            results = []

            # 进行dimer分析
            for i in range(len(sequences)):
                for j in range(i + 1, len(sequences)):
                    seq1 = dcP[sequences[i]]
                    seq2 = dcP[sequences[j]]
                    diTm = primer3.calc_heterodimer_tm(seq1, seq2)

                    if diTm > tmin:
                        dime = primer3.calc_heterodimer(seq1, seq2, output_structure=True)
                        results.append({
                            'seq1': sequences[i],
                            'seq2': sequences[j],
                            'Tm': diTm,
                            '结构': dime.ascii_structure
                        })

                # 按Tm值排序结果
            sorted_results = sorted(results, key=lambda x: x['Tm'], reverse=True)

                # 结果展示
            st.subheader("分析结果:")
            search_query = st.text_input("搜索结果:")

            for res in sorted_results:
                if search_query.lower() in res['seq1'].lower() or search_query.lower() in res['seq2'].lower():
                    st.write(f"seq1: {res['seq1']}, seq2: {res['seq2']}, Tm: {res['Tm']}")
                    st.write(f"结构:\n```\n{res['结构']}\n```")

# 右侧列：结果展示区（可以根据需求修改）
    st.empty()  # 可以在此添加更多内容
