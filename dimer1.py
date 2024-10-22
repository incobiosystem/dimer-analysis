import streamlit as st
import pandas as pd
import primer3
from Bio import Seq

# 文件上传
uploaded_file = st.file_uploader("选择Excel文件", type=["xlsx"])

# Tm 值输入
tmin = st.number_input("本次分析将输出高于以下温度的dimer (默认10℃):", value=10.0)

# 处理文件
if uploaded_file:
    df = pd.read_excel(uploaded_file)
    sequences = df['物料名称'].unique().tolist()  # 获取物料名称的唯一值
    selected_sequences = st.multiselect("选择序列进行分析:", sequences)

    # 开始分析按钮
    if st.button("开始分析"):
        if len(selected_sequences) < 2:
            st.warning("请至少选择两个序列进行分析！")
        else:
            dcP = {row['物料名称']: row['序列'] for _, row in df.iterrows() if row['物料名称'] in selected_sequences}

            results = []
            for i in range(len(selected_sequences)):
                for j in range(i + 1, len(selected_sequences)):
                    seq1 = dcP[selected_sequences[i]]
                    seq2 = dcP[selected_sequences[j]]

                    diTm = primer3.calc_heterodimer_tm(seq1, seq2)
                    if diTm > tmin:
                        dime = primer3.calc_heterodimer(seq1, seq2, output_structure=True)
                        results.append({
                            'seq1': selected_sequences[i],
                            'seq2': selected_sequences[j],
                            'Tm': diTm,
                            '结构': dime.ascii_structure
                        })

            # 按照Tm值排序
            sorted_results = sorted(results, key=lambda x: x['Tm'], reverse=True)

            # 显示结果
            if sorted_results:
                st.subheader("分析结果:")
                for res in sorted_results:
                    st.write(f"**seq1:** {res['seq1']}, **seq2:** {res['seq2']}, **Tm:** {res['Tm']}")
                    st.text(f"**结构:**\n{res['结构']}\n")
            else:
                st.write("没有符合条件的结果！")
