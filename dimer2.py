import streamlit as st
import pandas as pd
import primer3
from Bio import Seq
import os
import math
import itertools
import io

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

def process_sequences(df):
    """
    处理DataFrame中的序列，展开简并碱基。
    """
    # 存储输出数据
    output_data = []

    # 遍历每一行name和sequence
    for _, row in df.iterrows():
        if pd.isna(row['name']) or pd.isna(row['sequence']):
            continue  # 跳过不完整的数据
            
        material_name = row['name']
        sequence = row['sequence']
        
        # 检查序列是否包含简并碱基
        has_degenerate = any(base in degenerate_bases for base in sequence)
        
        if has_degenerate:
            # 展开所有可能的序列组合
            expanded_sequences = expand_degenerate_sequence(sequence)
            
            # 为每个生成的序列创建新的物料名称
            for i, seq in enumerate(expanded_sequences, start=1):
                new_material_name = f"{material_name}-{i}"
                output_data.append([new_material_name, seq])
        else:
            # 如果没有简并碱基，直接添加原始序列
            output_data.append([material_name, sequence])

    # 将结果输出为 DataFrame
    if output_data:
        output_df = pd.DataFrame(output_data, columns=['name', 'sequence'])
        return output_df
    else:
        return pd.DataFrame(columns=['name', 'sequence'])

# Matthews 2004 DNA能量模型参数
class ThermodynamicParams:
    # 悬空端(dangling ends)的热力学参数
    DANGLING_ENDS_DH = {
        'AT': -2.2,  # kcal/mol
        'TA': -2.3,
        'GC': -3.4,
        'CG': -3.2
    }
    
    DANGLING_ENDS_DS = {
        'AT': -6.9,  # cal/mol·K
        'TA': -7.2,
        'GC': -10.1,
        'CG': -9.8
    }

def adjust_thermodynamic_params(seq1, seq2, dg, tm, temp_c=37.0, dna_conc=50.0):
    """
    根据Matthews 2004模型调整热力学参数
    
    Args:
        seq1: 第一条序列
        seq2: 第二条序列
        dg: 原始deltaG值 (cal/mol)
        tm: 原始Tm值 (°C)
        temp_c: 模拟温度 (°C)
        
    Returns:
        adjusted_dg: 调整后的deltaG值 (cal/mol)
        adjusted_tm: 调整后的Tm值 (°C)
    """
    # 获取末端碱基对
    end5 = seq1[0] + seq2[-1]
    end3 = seq1[-1] + seq2[0]
    
    # 计算悬空端的贡献
    dh_dangling = ThermodynamicParams.DANGLING_ENDS_DH.get(end5, 0) + \
                 ThermodynamicParams.DANGLING_ENDS_DH.get(end3, 0)
    ds_dangling = ThermodynamicParams.DANGLING_ENDS_DS.get(end5, 0) + \
                 ThermodynamicParams.DANGLING_ENDS_DS.get(end3, 0)
    
    # 计算温度修正（所有计算都在cal/mol单位下进行）
    temp_k = temp_c + 273.15
    dh_dangling = dh_dangling * 1000  # 从kcal/mol转换到cal/mol
    dg_dangling = dh_dangling - (temp_k * ds_dangling)  # ds_dangling已经是cal/mol·K
    
    # 调整deltaG
    adjusted_dg = dg + dg_dangling
    
    # 调整Tm
    # 使用修正后的热力学参数计算Tm的变化
    # 使用完整的热力学方程：ΔTm = ΔH/(ΔS + R * ln(CT/4))
    R = 1.987  # 气体常数，cal/(mol·K)
    CT = dna_conc * 1e-9  # 将 nM 转换为 M
    
    if ds_dangling != 0:
        delta_tm = dh_dangling / (ds_dangling + R * math.log(CT/4))
        # 由于是温度变化，不需要加273.15
        adjusted_tm = tm + delta_tm * 0.1  # 添加一个缩放因子来减小影响
    else:
        adjusted_tm = tm
    
    return adjusted_dg, adjusted_tm

# 初始化session_state
if 'analysis_done' not in st.session_state:
    st.session_state.analysis_done = False
    st.session_state.sorted_results = []
    st.session_state.sequences = []
    st.session_state.expanded_df = None

# 设置页面标题和图标
st.set_page_config(page_title="Dimer Analysis with Degenerate Bases", page_icon="🔬", layout="wide")  # 使用宽布局

# 顶部布局，三列布局，带有间隔
col1, col_spacer1, col2, col_spacer2, col3 = st.columns([1.5, 0.5, 3, 0.5, 4.5])  # 使用间隔列来调整间距

# 左侧导航栏
with col1:
    logo_path = "logo.png"  # 确保logo文件名与路径正确
    if os.path.exists(logo_path):
        st.image(logo_path, width=300)  # 调整图片宽度
    st.markdown("<h1 style='font-size: 32px;'>Dimer Analysis</h1>", unsafe_allow_html=True)  # 可调整标题字体大小
    st.markdown("<h1 style='font-size: 16px;'>Incobiosystem©</h5>", unsafe_allow_html=True)  # 可调整标题字体大小
    
    # 添加计算逻辑简介在左下角
    st.markdown("---")
  
    st.markdown("""
    **应用原理**:
    1. 自动解析简并碱基并展开为所有可能的序列组合
    2. 使用primer3计算二聚体结构和初始热力学参数
    3. 应用Matthews 2004模型调整热力学参数
    4. 考虑悬空端(dangling ends)对热力学参数的贡献
    5. 计算调整后的Tm和ΔG值
    
    **参考文献**: Matthews et al. (2004) Biochemistry
    """)

# 中间列：上传文件和选择序列
with col2:
    st.subheader("上传需要分析的Excel文件")
    uploaded_file = st.file_uploader("选择要分析的文件", type=["xlsx"])
    # 输入 deltaG 阈值
    dg_threshold = st.number_input("分析将输出deltaG小于以下值的二聚体 (默认 0 cal/mol)", value=0.0)
    # 输入模拟温度
    simulation_temp = st.number_input("二聚体计算的模拟温度 (°C)", value=37.0, help="温度会影响二聚体结构的形成和deltaG值，但不直接影响Tm值的计算。Tm是由热力学参数计算得出的熔解温度，而deltaG是在指定温度下的自由能变化。")
    # 输入离子浓度参数
    with st.expander("高级参数设置", expanded=False):
        mv_conc = st.number_input("单价离子浓度 (mM)", value=50.0)
        dv_conc = st.number_input("二价离子浓度 (mM)", value=1.5)
        dntp_conc = st.number_input("dNTP浓度 (mM)", value=0.6)
        dna_conc = st.number_input("DNA浓度 (nM)", value=50.0)

    if uploaded_file is not None:
        df = pd.read_excel(uploaded_file)
        
        # 检查列名是否为name和sequence，如果不是，尝试使用物料名称和序列
        if 'name' not in df.columns or 'sequence' not in df.columns:
            if '物料名称' in df.columns and '序列' in df.columns:
                # 重命名列
                df = df.rename(columns={'物料名称': 'name', '序列': 'sequence'})
                st.info("已将列名'物料名称'和'序列'重命名为'name'和'sequence'")
            else:
                st.error("Excel文件必须包含'name'和'sequence'列或'物料名称'和'序列'列")
                st.stop()
        
        # 处理简并碱基
        if st.session_state.expanded_df is None:
            with st.spinner("正在处理简并碱基..."):
                expanded_df = process_sequences(df)
                st.session_state.expanded_df = expanded_df
                
                # 显示简并碱基展开后的序列数量
                original_count = len(df)
                expanded_count = len(expanded_df)
                if expanded_count > original_count:
                    st.success(f"已将{original_count}个含简并碱基的序列展开为{expanded_count}个序列")
                    
                    # 提供下载展开后的序列
                    buffer = io.BytesIO()
                    with pd.ExcelWriter(buffer, engine='openpyxl') as writer:
                        expanded_df.to_excel(writer, index=False)
                    
                    st.download_button(
                        label="下载展开后的序列文件",
                        data=buffer.getvalue(),
                        file_name=uploaded_file.name.replace('.xlsx', '_expanded.xlsx'),
                        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                    )
        else:
            expanded_df = st.session_state.expanded_df

        sequences = []
        st.subheader("选择要分析的序列:")
        # 添加全选按钮
        select_all = st.checkbox("全选/取消全选")
        
        # 创建一个字典来存储所有有效的物料名称
        valid_materials = {}
        for index, row in expanded_df.iterrows():
            if pd.isna(row['name']) or pd.isna(row['sequence']):
                st.warning(f"第 {index+1} 行数据不完整，跳过该行")
                continue  # 跳过不完整的数据
            valid_materials[row['name']] = row['sequence']
        
        # 根据全选按钮状态设置所有复选框
        for material_name in valid_materials.keys():
            selected = st.checkbox(material_name, value=select_all)
            if selected:
                sequences.append(material_name)

with col3:
    # 开始分析按钮
    if st.button("开始分析"):  # 可以调整按钮大小和样式
        if len(sequences) < 2:
            st.warning("请至少选择两个序列进行分析。")
        else:
            # 保存选择的序列到session_state
            st.session_state.sequences = sequences.copy()
            
            dcP = {row['name']: row['sequence'] for index, row in expanded_df.iterrows() if row['name'] in sequences}
            results = []

            # 进行dimer分析
            with st.spinner("正在进行二聚体分析..."):
                for i in range(len(sequences)):
                    seq1 = dcP[sequences[i]]
                    # 与自身序列做dimer分析
                    dime_self = primer3.calc_heterodimer(seq1, seq1, 
                                                        mv_conc=mv_conc, 
                                                        dv_conc=dv_conc, 
                                                        dntp_conc=dntp_conc, 
                                                        dna_conc=dna_conc, 
                                                        temp_c=simulation_temp,
                                                        output_structure=True)
                    if dime_self.dg < dg_threshold:
                        # 计算自身二聚体的Tm值
                        diTm_self = primer3.calc_heterodimer_tm(seq1, seq1, 
                                                              mv_conc=mv_conc, 
                                                              dv_conc=dv_conc, 
                                                              dntp_conc=dntp_conc, 
                                                              dna_conc=dna_conc)
                        
                        # 应用Matthews 2004模型调整热力学参数
                        adjusted_dg, adjusted_tm = adjust_thermodynamic_params(
                            seq1, seq1, dime_self.dg, diTm_self, simulation_temp, dna_conc)
                        
                        results.append({
                            'seq1': sequences[i],
                            'seq2': sequences[i],
                            'Tm': adjusted_tm,
                            'deltaG': adjusted_dg,
                            '结构': dime_self.ascii_structure
                        })

                    for j in range(i + 1, len(sequences)):
                        seq2 = dcP[sequences[j]]
                        # 计算两序列间二聚体的Tm值，不传入simulation_temp参数
                        diTm = primer3.calc_heterodimer_tm(seq1, seq2, 
                                                         mv_conc=mv_conc, 
                                                         dv_conc=dv_conc, 
                                                         dntp_conc=dntp_conc, 
                                                         dna_conc=dna_conc)

                        dimer = primer3.calc_heterodimer(seq1, seq2, 
                                                          mv_conc=mv_conc, 
                                                          dv_conc=dv_conc, 
                                                          dntp_conc=dntp_conc, 
                                                          dna_conc=dna_conc, 
                                                          temp_c=simulation_temp,
                                                          output_structure=True)
                        if dimer.dg < dg_threshold:
                            # 应用Matthews 2004模型调整热力学参数
                            adjusted_dg, adjusted_tm = adjust_thermodynamic_params(
                                seq1, seq2, dimer.dg, diTm, simulation_temp, dna_conc)
                            
                            results.append({
                                'seq1': sequences[i],
                                'seq2': sequences[j],
                                'Tm': adjusted_tm,
                                'deltaG': adjusted_dg,
                                '结构': dimer.ascii_structure
                            })

            # 按deltaG值排序结果（由小到大）
            sorted_results = sorted(results, key=lambda x: x['deltaG'])
            
            # 保存分析结果到session_state
            st.session_state.sorted_results = sorted_results
            st.session_state.analysis_done = True
            
            st.success(f"分析完成！找到 {len(sorted_results)} 个二聚体结构。")
    
    # 如果已经完成分析，显示结果
    if st.session_state.analysis_done:
        # 结果展示
        st.subheader("分析结果:")
        col1, col2 = st.columns(2)
        with col1:
            search_query = st.text_input("搜索结果:")
        with col2:
            # 添加序列筛选下拉框
            filter_sequence = st.selectbox(
                "筛选特定序列的二聚体:",
                ["所有序列"] + st.session_state.sequences,
                key="filter_sequence_selector"
            )

        # 根据搜索条件和序列筛选条件筛选结果
        filtered_results = []
        for res in st.session_state.sorted_results:
            # 首先检查序列筛选条件
            if filter_sequence == "所有序列" or filter_sequence == res['seq1'] or filter_sequence == res['seq2']:
                # 然后检查搜索条件
                if not search_query or search_query.lower() in res['seq1'].lower() or search_query.lower() in res['seq2'].lower():
                    filtered_results.append(res)
        
        # 显示筛选结果数量
        st.write(f"找到 {len(filtered_results)} 个符合条件的二聚体结构")
        
        # 显示筛选后的结果
        for res in filtered_results:
            with st.container():
                # 使用更紧凑的布局
                st.markdown(f"**序列1**: {res['seq1']} | **序列2**: {res['seq2']}")
                col1, col2 = st.columns([3, 7])
                with col1:
                    st.markdown(f"Tm={res['Tm']:.2f}°C, ΔG={res['deltaG']:.2f} cal/mol")
                with col2:
                    st.code(res['结构'], language=None)
                st.markdown("---")

    # 右侧列：结果展示区（可以根据需求修改）
    st.empty()  # 可以在此添加更多内容
