import pandas as pd
import itertools
import streamlit as st

# 简并碱基与其对应的可能碱基的映射
degenerate_bases = {
    'R': ['A', 'G'],
    'Y': ['C', 'T'],
    'S': ['G', 'C'],
    'W': ['A', 'T'],
    'K': ['G', 'T'],
    'M': ['A', 'C'],
    'B': ['C', 'G', 'T'],
    'D': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'],
    'V': ['A', 'C', 'G'],
    'N': ['A', 'C', 'G', 'T']
}

def expand_degenerate_sequence(sequence):
    """
    展开包含简并碱基的引物序列，返回所有可能的组合。
    """
    # 生成一个列表，包含序列中每个碱基对应的替换选项
    replacement_options = []
    for base in sequence:
        if base in degenerate_bases:
            replacement_options.append(degenerate_bases[base])
        else:
            replacement_options.append([base])  # 普通碱基，直接保留

    # 使用 itertools.product 生成所有组合
    expanded_sequences = [''.join(comb) for comb in itertools.product(*replacement_options)]
    return expanded_sequences

def analyze_degenerate_sequence(sequence):
    """
    分析序列中简并碱基的数量，以及每个简并碱基代表的可能碱基种类。
    """
    degenerate_count = 0
    degenerate_details = {}

    for base in sequence:
        if base in degenerate_bases:
            degenerate_count += 1
            degenerate_details[base] = degenerate_bases[base]

    return degenerate_count, degenerate_details

def process_sequences(uploaded_file):
    """
    处理上传的 Excel 文件，输出处理后的结果。
    """
    df = pd.read_excel(uploaded_file)

    # 存储输出数据
    output_data = []

    # 遍历每一行物料名称和序列
    for _, row in df.iterrows():
        material_name = row['物料名称']
        sequence = row['序列']

        # 展开所有可能的序列组合
        expanded_sequences = expand_degenerate_sequence(sequence)

        # 为每个生成的序列创建新的物料名称
        for i, seq in enumerate(expanded_sequences, start=1):
            new_material_name = f"{material_name}-{i}"
            output_data.append([new_material_name, seq])

    # 将结果输出为 DataFrame
    output_df = pd.DataFrame(output_data, columns=['物料名称', '序列'])
    return output_df

def main():
    st.title('引物序列简并碱基转化工具')

    # 文件上传控件
    uploaded_file = st.file_uploader("上传包含序列的Excel文件", type="xlsx")

    if uploaded_file is not None:
        # 处理上传的文件
        output_df = process_sequences(uploaded_file)
        
        # 显示处理后的 DataFrame
        st.write(output_df)

        # 允许用户下载处理后的 Excel 文件
        output_excel = 'output_sequences.xlsx'
        output_df.to_excel(output_excel, index=False)
        
        st.download_button(
            label="下载处理后的文件",
            data=open(output_excel, "rb").read(),
            file_name=output_excel,
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        )

if __name__ == "__main__":
    main()
