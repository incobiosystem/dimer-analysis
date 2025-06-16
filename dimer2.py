import streamlit as st
import pandas as pd
import primer3
from Bio import Seq
import os
import math
import itertools
import io
from itertools import combinations
import streamlit as st
from datetime import datetime

def calculate_combination_score(sequences, target_map, dimer_results, max_acceptable_deltaG=-500.0, min_acceptable_deltaG=-5000.0):
    """
    计算一个组合中所有二聚体的deltaG统计信息
    
    Args:
        sequences: 序列名称列表
        target_map: 序列名称到靶标的映射
        dimer_results: 所有二聚体分析结果
        max_acceptable_deltaG: 可接受的最大deltaG值（最负值）
        min_acceptable_deltaG: 可接受的最小deltaG值（最负值）
    
    Returns:
        dict: 包含deltaG统计信息的字典
    """
    combination_dimers = []
    problematic_dimers = []
    
    # 找出这个组合中所有相关的二聚体
    for result in dimer_results:
        if result['seq1'] in sequences and result['seq2'] in sequences:
            combination_dimers.append(result['deltaG'])
            # 检查是否低于用户设定的最小可接受deltaG值
            if result['deltaG'] < min_acceptable_deltaG:
                problematic_dimers.append({
                    'seq1': result['seq1'],
                    'seq2': result['seq2'],
                    'deltaG': result['deltaG'],
                    'problem_type': '不符合要求'
                })
    
    if not combination_dimers:
        return {
            'count': 0,
            'total': 0,
            'min': 0,
            'max': 0,
            'average': 0,
            'sequences': sequences,
            'targets': [target_map.get(seq, '未知') for seq in sequences],
            'problematic_dimers': [],
            'problematic_count': 0
        }
    
    return {
        'count': len(combination_dimers),
        'total': sum(combination_dimers),
        'min': min(combination_dimers),
        'max': max(combination_dimers),
        'average': sum(combination_dimers) / len(combination_dimers),
        'sequences': sequences,
        'targets': [target_map.get(seq, '未知') for seq in sequences],
        'problematic_dimers': problematic_dimers,
        'problematic_count': len(problematic_dimers)
    }

def generate_target_combinations(target_sequences, min_combination_size=1, max_combination_size=5):
    """
    生成靶标组合建议
    
    Args:
        target_sequences: 靶标到序列列表的映射
        min_combination_size: 最小组合大小
        max_combination_size: 最大组合大小
    
    Returns:
        list: 所有可能的靶标组合
    """
    targets = list(target_sequences.keys())
    all_combinations = []
    
    # 生成min_combination_size到max_combination_size大小的所有组合
    for size in range(min_combination_size, min(len(targets) + 1, max_combination_size + 1)):
        for combo in combinations(targets, size):
            # 获取这个靶标组合中的所有序列
            sequences = []
            for target in combo:
                sequences.extend(target_sequences[target])
            
            all_combinations.append({
                'targets': list(combo),
                'sequences': sequences
            })
    
    return all_combinations

def intelligent_well_assignment(targets, target_sequences, dimer_results, target_map, min_acceptable_deltaG=-3000.0, min_targets_per_well=1, max_targets_per_well=4, max_wells=4, max_iterations=1000000):
    """
    智能分配靶标到反应孔的算法
    
    Args:
        targets: 所有靶标列表
        target_sequences: 靶标到序列的映射
        dimer_results: 二聚体分析结果
        target_map: 序列到靶标的映射
        max_acceptable_deltaG: 可接受的最大deltaG值（最负值）
        min_targets_per_well: 单孔最小靶标数量
        max_targets_per_well: 单孔最大靶标数量
        max_wells: 最大孔数
    
    Returns:
        list: 最优的孔分配方案
    """
    import random
    from itertools import combinations
    
    # 创建dimer查找字典，加速查找
    dimer_lookup = {}
    for result in dimer_results:
        key1 = (result['seq1'], result['seq2'])
        key2 = (result['seq2'], result['seq1'])
        dimer_lookup[key1] = result['deltaG']
        dimer_lookup[key2] = result['deltaG']
    
    def get_dimer_deltaG(seq1, seq2):
        """获取两个序列间的dimer deltaG值"""
        return dimer_lookup.get((seq1, seq2), 0)
    
    def evaluate_well(target_list):
        """评估一个孔内靶标组合的质量"""
        sequences = []
        for target in target_list:
            sequences.extend(target_sequences[target])
        
        worst_deltaG = 0
        problematic_count = 0
        all_deltaGs = []
        
        # 计算所有序列对的dimer deltaG
        for i in range(len(sequences)):
            for j in range(i, len(sequences)):
                deltaG = get_dimer_deltaG(sequences[i], sequences[j])
                all_deltaGs.append(deltaG)
                
                if deltaG < min_acceptable_deltaG:
                    problematic_count += 1
                    if deltaG < worst_deltaG:
                        worst_deltaG = deltaG
        
        avg_deltaG = sum(all_deltaGs) / len(all_deltaGs) if all_deltaGs else 0
        
        return {
            'targets': target_list,
            'sequences': sequences,
            'worst_deltaG': worst_deltaG,
            'avg_deltaG': avg_deltaG,
            'problematic_count': problematic_count,
            'is_valid': problematic_count == 0,
            'score': problematic_count * 1000 + abs(worst_deltaG)  # 惩罚函数
        }
    
    def generate_well_assignments(targets, num_wells, max_combinations=200000):
        """生成将靶标分配到指定数量孔的所有可能方案（优化版：基于deltaG质量分配）"""
        if num_wells == 1:
            # 允许最后一个孔包含任意数量的靶标（包括少于最小值的情况）
            # 这样可以支持不均匀分配，如9个靶标分成4-4-1
            if len(targets) <= max_targets_per_well:
                return [[targets]]
            else:
                return []
        
        assignments = []
        # 尝试不同的第一个孔的靶标组合（从min_targets_per_well个到max_targets_per_well个）
        for well_size in range(min_targets_per_well, min(max_targets_per_well + 1, len(targets) - num_wells + 2)):
            # 添加剪枝：如果当前组合数已经足够，停止生成更多组合
            if len(assignments) >= max_combinations:
                break
                
            for first_well_targets in combinations(targets, well_size):
                remaining_targets = [t for t in targets if t not in first_well_targets]
                
                # 确保剩余靶标数量不超过剩余孔数的最大容量
                if len(remaining_targets) <= (num_wells - 1) * max_targets_per_well:
                    # 递归分配剩余靶标
                    sub_assignments = generate_well_assignments(remaining_targets, num_wells - 1, max_combinations - len(assignments))
                    
                    for sub_assignment in sub_assignments:
                        full_assignment = [list(first_well_targets)] + sub_assignment
                        assignments.append(full_assignment)
                        
                        # 添加剪枝：如果已经生成足够的组合，停止
                        if len(assignments) >= max_combinations:
                            break
                    
                    # 如果已经生成足够的组合，跳出循环
                    if len(assignments) >= max_combinations:
                        break
        
        return assignments
    
    # 生成所有可能的分配方案（优化版：基于deltaG质量评估）
    best_assignments = []
    calculation_count = 0
    
    # 添加性能优化：限制计算量
    max_assignments_per_well_count = 200000  # 每个孔数最多评估200000个方案
    
    for num_wells in range(1, min(max_wells + 1, len(targets) + 1)):
        assignments = generate_well_assignments(targets, num_wells)
        
        # 限制每个孔数的方案数量，避免组合爆炸
        if len(assignments) > max_assignments_per_well_count:
            import random
            assignments = random.sample(assignments, max_assignments_per_well_count)
        
        for assignment in assignments:
            calculation_count += 1
            
            # 早期终止条件：如果计算量超过限制，停止计算
            if calculation_count > max_iterations:
                break
                
            # 评估每个孔
            well_evaluations = []
            total_avg_deltaG = 0
            worst_deltaG_overall = 0
            all_valid = True
            total_problematic = 0
            
            for well_targets in assignment:
                well_eval = evaluate_well(well_targets)
                well_evaluations.append(well_eval)
                total_avg_deltaG += abs(well_eval['avg_deltaG'])
                if well_eval['worst_deltaG'] < worst_deltaG_overall:
                    worst_deltaG_overall = well_eval['worst_deltaG']
                if not well_eval['is_valid']:
                    all_valid = False
                total_problematic += well_eval['problematic_count']
            
            # 新的评分策略：优先考虑deltaG质量
            avg_deltaG_score = total_avg_deltaG / len(well_evaluations) if well_evaluations else 0
            worst_deltaG_penalty = abs(worst_deltaG_overall) * 0.1
            quality_score = avg_deltaG_score + worst_deltaG_penalty + total_problematic * 100
            
            assignment_result = {
                'wells': well_evaluations,
                'num_wells': num_wells,
                'total_score': quality_score,  # 新的质量评分
                'avg_deltaG_score': avg_deltaG_score,
                'worst_deltaG_overall': worst_deltaG_overall,
                'all_valid': all_valid,
                'total_problematic': total_problematic
            }
            
            best_assignments.append(assignment_result)
        
        # 如果已经达到计算限制，跳出外层循环
        if calculation_count > max_iterations:
            break
    
    # 存储计算次数到streamlit session state
    try:
        import streamlit as st
        st.session_state.total_calculations = calculation_count
    except:
        pass
    
    # 调试信息：打印生成的分配方案数量
    try:
        import streamlit as st
        st.write(f"🔍 调试信息：生成了 {len(best_assignments)} 个分配方案")
        if len(best_assignments) == 0:
            st.write(f"⚠️ 未生成任何分配方案，可能的原因：")
            st.write(f"   - 靶标数量: {len(targets)}")
            st.write(f"   - 最小靶标/孔: {min_targets_per_well}")
            st.write(f"   - 最大靶标/孔: {max_targets_per_well}")
            st.write(f"   - 最大孔数: {max_wells}")
            st.write(f"   - 理论最少需要孔数: {math.ceil(len(targets) / max_targets_per_well)}")
            if max_wells < math.ceil(len(targets) / max_targets_per_well):
                st.error(f"⚠️ 最大孔数({max_wells})小于理论最少需要孔数({math.ceil(len(targets) / max_targets_per_well)})，请增加最大孔数！")
    except:
        pass
    
    # 排序：优先选择所有孔都有效的方案，然后按deltaG质量排序
    best_assignments.sort(key=lambda x: (
        not x['all_valid'],  # 优先选择所有孔都有效的方案
        x['total_problematic'],  # 其次选择问题最少的方案
        -x['avg_deltaG_score'],  # 修正为数值最大（非绝对值）
        -x['worst_deltaG_overall']  # 修正为数值最大（非绝对值）
    ))
    
    return best_assignments[:10]  # 返回前10个最优方案

def optimize_combinations_with_ortools(expanded_df, dimer_results, min_acceptable_deltaG=-5000.0, min_targets_per_well=1, max_targets_per_well=6, design_wells=8):
    """
    使用Google OR-Tools的CP-SAT求解器优化靶标组合
    
    Args:
        expanded_df: 展开后的序列DataFrame
        dimer_results: 二聚体分析结果
        min_acceptable_deltaG: 最小可接受deltaG值（用于判断问题dimer）
        min_targets_per_well: 单孔最小靶标数量
        max_targets_per_well: 单孔最大靶标数量
        design_wells: 设计孔数（严格按此孔数分组）
    
    Returns:
        list: 最优的组合建议
    """
    try:
        from ortools.sat.python import cp_model
        import streamlit as st
        import time
    except ImportError:
        st.error("请安装Google OR-Tools: pip install ortools")
        return []
    
    if 'target' not in expanded_df.columns:
        return []
    
    start_time = time.time()
    
    # 只处理勾选的序列
    if hasattr(st.session_state, 'sequences') and st.session_state.sequences:
        # 过滤出勾选序列对应的数据
        filtered_df = expanded_df[expanded_df['name'].isin(st.session_state.sequences)]
    else:
        # 如果没有勾选序列，使用所有数据
        filtered_df = expanded_df
    
    # 创建靶标到序列的映射
    target_sequences = {}
    target_map = {}
    
    for _, row in filtered_df.iterrows():
        target = row.get('target', '未知')
        sequence = row['name']
        target_map[sequence] = target
        
        if target not in target_sequences:
            target_sequences[target] = []
        target_sequences[target].append(sequence)
    
    targets = list(target_sequences.keys())
    num_targets = len(targets)
    
    if num_targets == 0:
        st.warning("没有找到勾选的序列对应的靶标数据")
        return []
    
    # 创建dimer查找字典
    dimer_lookup = {}
    for result in dimer_results:
        key1 = (result['seq1'], result['seq2'])
        key2 = (result['seq2'], result['seq1'])
        dimer_lookup[key1] = result['deltaG']
        dimer_lookup[key2] = result['deltaG']
    
    def get_dimer_deltaG(seq1, seq2):
        """获取两个序列间的dimer deltaG值"""
        return dimer_lookup.get((seq1, seq2), 0)
    
    def calculate_well_problematic_deltaG_sum(target_indices):
        """计算一个孔内所有问题dimer的deltaG总和"""
        sequences = []
        for target_idx in target_indices:
            sequences.extend(target_sequences[targets[target_idx]])
        
        problematic_deltaG_sum = 0
        for i in range(len(sequences)):
            for j in range(i, len(sequences)):
                deltaG = get_dimer_deltaG(sequences[i], sequences[j])
                if deltaG < min_acceptable_deltaG:  # 问题dimer
                    problematic_deltaG_sum += deltaG
        
        return problematic_deltaG_sum
    
    # 预计算所有可能的孔组合的问题dimer deltaG总和
    st.write("🔍 正在预计算所有可能的孔组合...")
    well_combinations = {}
    calculation_count = 0
    size_counts = {}  # 记录每种大小的组合数量
    
    from itertools import combinations
    for well_size in range(min_targets_per_well, max_targets_per_well + 1):
        size_count = 0
        for target_combo in combinations(range(num_targets), well_size):
            deltaG_sum = calculate_well_problematic_deltaG_sum(target_combo)
            well_combinations[target_combo] = deltaG_sum
            calculation_count += 1
            size_count += 1
        size_counts[well_size] = size_count
    
    # 显示计算统计信息
    size_info = ", ".join([f"{size}个靶标: {count}种组合" for size, count in size_counts.items()])
    st.write(f"📊 预计算完成：总共 {calculation_count} 种组合 ({size_info})")
    
    # 创建CP-SAT模型
    model = cp_model.CpModel()
    
    # 决策变量：x[i][j] = 1 表示靶标i分配到孔j
    x = {}
    for i in range(num_targets):
        for j in range(design_wells):
            x[i, j] = model.NewBoolVar(f'x_{i}_{j}')
    
    # 约束1：每个靶标必须分配到且仅分配到一个孔
    for i in range(num_targets):
        model.Add(sum(x[i, j] for j in range(design_wells)) == 1)
    
    # 约束2：每个孔的靶标数量限制
    for j in range(design_wells):
        model.Add(sum(x[i, j] for i in range(num_targets)) >= min_targets_per_well)
        model.Add(sum(x[i, j] for i in range(num_targets)) <= max_targets_per_well)
    
    # 目标函数：最大化所有孔的问题dimer deltaG总和
    # 由于CP-SAT只支持整数，我们将deltaG值乘以1000并取整
    objective_terms = []
    
    # 为每个孔和每个可能的靶标组合创建指示变量
    for j in range(design_wells):
        for target_combo, deltaG_sum in well_combinations.items():
            combo_size = len(target_combo)
            if min_targets_per_well <= combo_size <= max_targets_per_well:
                # 创建指示变量：当且仅当孔j包含exactly这个靶标组合时为1
                combo_indicator = model.NewBoolVar(f'combo_{j}_{"_".join(map(str, sorted(target_combo)))}')
                
                # 约束：如果combo_indicator为1，则孔j必须包含target_combo中的所有靶标且不包含其他靶标
                # 孔j包含target_combo中的所有靶标
                for target_idx in target_combo:
                    model.Add(x[target_idx, j] >= combo_indicator)
                
                # 孔j不包含target_combo之外的靶标
                for target_idx in range(num_targets):
                    if target_idx not in target_combo:
                        model.Add(x[target_idx, j] <= 1 - combo_indicator)
                
                # 如果孔j恰好包含target_combo中的所有靶标且不包含其他靶标，则combo_indicator为1
                model.Add(sum(x[target_idx, j] for target_idx in target_combo) >= combo_size * combo_indicator)
                model.Add(sum(x[target_idx, j] for target_idx in range(num_targets)) <= combo_size + (num_targets - combo_size) * (1 - combo_indicator))
                
                # 将deltaG贡献添加到目标函数
                deltaG_contribution = int(deltaG_sum * 1000)  # 转换为整数
                objective_terms.append(combo_indicator * deltaG_contribution)
    
    # 设置目标函数：最大化问题dimer deltaG总和（找到最接近0的值，即绝对值最小）
    model.Maximize(sum(objective_terms))
    
    # 求解器配置：确保确定性结果
    solver = cp_model.CpSolver()
    solver.parameters.random_seed = 42  # 设置固定随机种子
    solver.parameters.num_search_workers = 1  # 禁用并行化
    solver.parameters.max_time_in_seconds = 120.0  # 增加求解时间
    solver.parameters.log_search_progress = False  # 禁用日志输出
    
    st.write("🚀 正在使用OR-Tools求解全局最优方案...")
    status = solver.Solve(model)
    
    end_time = time.time()
    calculation_time = end_time - start_time
    
    # 保存计算统计信息到session_state
    st.session_state.calculation_time = calculation_time
    st.session_state.calculation_count = calculation_count
    
    if status == cp_model.OPTIMAL:
        objective_value = solver.ObjectiveValue() / 1000
        st.success(f"✅ 找到全局最优解！目标函数值（总体问题deltaG，越接近0越好）: {objective_value:.1f} cal/mol")
        
        # 提取解决方案
        solution_wells = [[] for _ in range(design_wells)]
        for i in range(num_targets):
            for j in range(design_wells):
                if solver.Value(x[i, j]) == 1:
                    solution_wells[j].append(targets[i])
        
        # 计算实际的总体deltaG（验证）
        total_problematic_deltaG = 0
        for well_targets in solution_wells:
            if well_targets:
                sequences = []
                for target in well_targets:
                    sequences.extend(target_sequences[target])
                
                for i in range(len(sequences)):
                    for j in range(i, len(sequences)):
                        deltaG = get_dimer_deltaG(sequences[i], sequences[j])
                        if deltaG < min_acceptable_deltaG:
                            total_problematic_deltaG += deltaG
        
        st.info(f"📊 验证：实际总体问题deltaG = {total_problematic_deltaG:.1f} cal/mol")
        
        # 转换为原有格式
        scored_combinations = []
        for j, well_targets in enumerate(solution_wells):
            if well_targets:  # 只处理非空孔
                sequences = []
                for target in well_targets:
                    sequences.extend(target_sequences[target])
                
                score = calculate_combination_score(
                    sequences, 
                    target_map, 
                    dimer_results,
                    -3000.0,  # max_acceptable_deltaG固定值
                    min_acceptable_deltaG
                )
                score['target_names'] = well_targets
                score['well_number'] = j + 1
                score['assignment_id'] = 0
                score['total_wells'] = design_wells
                score['assignment_valid'] = score['problematic_count'] == 0
                score['total_problematic_deltaG'] = total_problematic_deltaG  # 添加总体deltaG
                scored_combinations.append(score)
        
        return scored_combinations
    
    elif status == cp_model.FEASIBLE:
        objective_value = solver.ObjectiveValue() / 1000
        st.warning(f"⚠️ 找到可行解（未达到全局最优）！目标函数值（总体问题deltaG，越接近0越好）: {objective_value:.1f} cal/mol")
        
        # 提取解决方案（与最优解处理相同）
        solution_wells = [[] for _ in range(design_wells)]
        for i in range(num_targets):
            for j in range(design_wells):
                if solver.Value(x[i, j]) == 1:
                    solution_wells[j].append(targets[i])
        
        # 计算实际的总体deltaG
        total_problematic_deltaG = 0
        for well_targets in solution_wells:
            if well_targets:
                sequences = []
                for target in well_targets:
                    sequences.extend(target_sequences[target])
                
                for i in range(len(sequences)):
                    for j in range(i, len(sequences)):
                        deltaG = get_dimer_deltaG(sequences[i], sequences[j])
                        if deltaG < min_acceptable_deltaG:
                            total_problematic_deltaG += deltaG
        
        st.info(f"📊 验证：实际总体问题deltaG = {total_problematic_deltaG:.1f} cal/mol")
        
        scored_combinations = []
        for j, well_targets in enumerate(solution_wells):
            if well_targets:
                sequences = []
                for target in well_targets:
                    sequences.extend(target_sequences[target])
                
                score = calculate_combination_score(
                    sequences, 
                    target_map, 
                    dimer_results,
                    -3000.0,
                    min_acceptable_deltaG
                )
                score['target_names'] = well_targets
                score['well_number'] = j + 1
                score['assignment_id'] = 0
                score['total_wells'] = design_wells
                score['assignment_valid'] = score['problematic_count'] == 0
                score['total_problematic_deltaG'] = total_problematic_deltaG  # 添加总体deltaG
                scored_combinations.append(score)
        
        return scored_combinations
    
    else:
        st.error("❌ 未找到可行解，请调整参数设置（如增加孔数或调整靶标数量限制）")
        return []

def optimize_combinations(expanded_df, dimer_results, max_combinations=10, max_acceptable_deltaG=-500.0, min_acceptable_deltaG=-5000.0, min_targets_per_well=1, max_targets_per_well=6, max_wells=8):
    """
    使用OR-Tools优化靶标组合（仅最优算法）
    
    Args:
        expanded_df: 展开后的序列DataFrame
        dimer_results: 二聚体分析结果
        max_combinations: 返回的最大组合数量（保留兼容性，实际不使用）
        max_acceptable_deltaG: 可接受的最大deltaG值（保留兼容性，实际不使用）
        min_acceptable_deltaG: 可接受的最小deltaG值（最负值）
        min_targets_per_well: 单孔最小靶标数量
        max_targets_per_well: 单孔最大靶标数量
        max_wells: 设计孔数（严格按此孔数分组）
    
    Returns:
        list: 全局最优的组合建议
    """
    # 直接使用OR-Tools求解，不再回退到传统算法
    return optimize_combinations_with_ortools(
        expanded_df, 
        dimer_results, 
        min_acceptable_deltaG, 
        min_targets_per_well, 
        max_targets_per_well, 
        max_wells
    )

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
    
    # 检查是否有靶标列
    has_target = 'target' in df.columns

    # 遍历每一行name和sequence
    for _, row in df.iterrows():
        if pd.isna(row['name']) or pd.isna(row['sequence']):
            continue  # 跳过不完整的数据
            
        material_name = row['name']
        sequence = row['sequence']
        target = row.get('target', '') if has_target else ''
        
        # 检查序列是否包含简并碱基
        has_degenerate = any(base in degenerate_bases for base in sequence)
        
        if has_degenerate:
            # 展开所有可能的序列组合
            expanded_sequences = expand_degenerate_sequence(sequence)
            
            # 为每个生成的序列创建新的物料名称
            for i, seq in enumerate(expanded_sequences, start=1):
                new_material_name = f"{material_name}-{i}"
                if has_target:
                    output_data.append([new_material_name, seq, target])
                else:
                    output_data.append([new_material_name, seq])
        else:
            # 如果没有简并碱基，直接添加原始序列
            if has_target:
                output_data.append([material_name, sequence, target])
            else:
                output_data.append([material_name, sequence])

    # 将结果输出为 DataFrame
    if output_data:
        if has_target:
            output_df = pd.DataFrame(output_data, columns=['name', 'sequence', 'target'])
        else:
            output_df = pd.DataFrame(output_data, columns=['name', 'sequence'])
        return output_df
    else:
        if has_target:
            return pd.DataFrame(columns=['name', 'sequence', 'target'])
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
    
    # 计算逻辑简介
    st.markdown("**🔬 计算逻辑简介**")
    st.markdown("""
    1. **简并碱基展开**: 自动识别并展开IUPAC简并碱基
    2. **二聚体分析**: 使用Primer3算法计算所有序列对的二聚体结构
    3. **智能合孔**: 基于ΔG阈值优化靶标分配到反应孔
    """)
  


# 初始化sequences变量在更外层作用域
sequences = []

# 中间列：上传文件和选择序列
with col2:
    st.markdown("<h3 style='font-size: 18px;'>📁 上传需要分析的Excel或CSV文件</h3>", unsafe_allow_html=True)
    uploaded_file = st.file_uploader("选择要分析的文件", type=["xlsx", "csv"])
    # 输入 deltaG 阈值
    dg_threshold = st.number_input("分析将输出deltaG小于以下值的二聚体 (默认 -3000 cal/mol)", value=-3000.0)
    # 输入模拟温度
    simulation_temp = st.number_input("二聚体计算的模拟温度 (°C)", value=37.0, help="温度会影响二聚体结构的形成和deltaG值，但不直接影响Tm值的计算。Tm是由热力学参数计算得出的熔解温度，而deltaG是在指定温度下的自由能变化。")
    
    # 输入离子浓度参数
    with st.expander("高级参数设置", expanded=False):
        mv_conc = st.number_input("单价离子浓度 (mM)", value=50.0)
        dv_conc = st.number_input("二价离子浓度 (mM)", value=1.5)
        dntp_conc = st.number_input("dNTP浓度 (mM)", value=0.6)
        dna_conc = st.number_input("DNA浓度 (nM)", value=50.0)

    if uploaded_file is not None:
        # 根据文件类型读取数据
        file_extension = uploaded_file.name.split('.')[-1].lower()
        
        if file_extension == 'csv':
            df = pd.read_csv(uploaded_file)
        else:  # xlsx
            df = pd.read_excel(uploaded_file)
        
        # 将列名转换为小写以便不区分大小写比较
        df.columns = [col.lower() for col in df.columns]
        
        # 创建列名映射字典
        name_columns = ['name', 'material name', 'materialname', '物料名称', '物料', '名称']
        sequence_columns = ['sequence', 'seq', '序列']
        target_columns = ['target', 'pathogen', '靶标', '病原体', '检测目标', '目标']
        
        # 检查是否存在name和sequence列（不区分大小写）
        name_col = None
        for col in name_columns:
            if col in df.columns:
                name_col = col
                break
                
        seq_col = None
        for col in sequence_columns:
            if col in df.columns:
                seq_col = col
                break
        
        # 检查是否存在靶标列
        target_col = None
        for col in target_columns:
            if col in df.columns:
                target_col = col
                break
        
        if name_col is not None and seq_col is not None:
            # 如果列名不是标准的'name'和'sequence'，则重命名
            rename_dict = {}
            if name_col != 'name':
                rename_dict[name_col] = 'name'
            if seq_col != 'sequence':
                rename_dict[seq_col] = 'sequence'
            if target_col is not None and target_col != 'target':
                rename_dict[target_col] = 'target'
            
            if rename_dict:
                df = df.rename(columns=rename_dict)
                renamed_cols = ', '.join([f"'{old}'->''{new}'" for old, new in rename_dict.items()])
                st.info(f"已重命名列: {renamed_cols}")
        else:
            st.error("文件必须包含name/物料名称和sequence/序列列（不区分大小写）")
            st.stop()
        
        # 检查是否有靶标信息
        has_target_info = 'target' in df.columns and not df['target'].isna().all()
        
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
                    
                    # 根据原始文件类型决定下载文件的格式
                    if file_extension == 'csv':
                        csv_data = expanded_df.to_csv(index=False)
                        buffer.write(csv_data.encode())
                        download_filename = uploaded_file.name.replace('.csv', '_expanded.csv')
                        mime_type = "text/csv"
                    else:  # xlsx
                        with pd.ExcelWriter(buffer, engine='openpyxl') as writer:
                            expanded_df.to_excel(writer, index=False)
                        download_filename = uploaded_file.name.replace('.xlsx', '_expanded.xlsx')
                        mime_type = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                    
                    buffer.seek(0)
                    
                    st.download_button(
                        label="下载展开后的序列文件",
                        data=buffer.getvalue(),
                        file_name=download_filename,
                        mime=mime_type
                    )
        else:
            expanded_df = st.session_state.expanded_df

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
        # 检查是否有靶标信息
        has_target_info = 'target' in st.session_state.expanded_df.columns and not st.session_state.expanded_df['target'].isna().all()
        
        # 如果有靶标信息，显示智能合孔方案
        if has_target_info:
            st.markdown("<h3 style='font-size: 18px;'>🎯 智能合孔方案</h3>", unsafe_allow_html=True)
            
            # 获取勾选序列对应的靶标
            if hasattr(st.session_state, 'sequences') and st.session_state.sequences:
                # 过滤出勾选序列对应的数据
                selected_df = st.session_state.expanded_df[st.session_state.expanded_df['name'].isin(st.session_state.sequences)]
                selected_targets = selected_df['target'].unique().tolist()
                st.info(f"检测到 {len(selected_targets)} 个靶标")
            else:
                # 如果没有勾选序列，使用所有靶标
                all_targets = st.session_state.expanded_df['target'].unique().tolist()
                st.info(f"检测到 {len(all_targets)} 个靶标")
            
            # 组合设计参数
            st.markdown("**⚙️ 组合设计参数**")
            col_param1, col_param2 = st.columns(2)
            
            with col_param1:
                min_acceptable_deltaG = st.number_input(
                    "最小可接受ΔG值 (cal/mol)", 
                    min_value=-50000.0, 
                    max_value=-100.0, 
                    value=-5000.0, 
                    step=100.0,
                    help="组合中二聚体ΔG的最小可接受值（最负值）"
                )
                
                max_targets_per_well = st.number_input(
                    "单孔最大靶标数量", 
                    min_value=1, 
                    max_value=20, 
                    value=4, 
                    step=1,
                    help="每个反应孔中允许的最大靶标数量"
                )
            
            with col_param2:
                min_targets_per_well = st.number_input(
                    "单孔最小靶标数量", 
                    min_value=1, 
                    max_value=10, 
                    value=3, 
                    step=1,
                    help="每个反应孔中要求的最小靶标数量"
                )
                
                # 获取当前可用的靶标数量来设置默认值
                if hasattr(st.session_state, 'sequences') and st.session_state.sequences:
                    selected_df = st.session_state.expanded_df[st.session_state.expanded_df['name'].isin(st.session_state.sequences)]
                    current_targets = selected_df['target'].unique().tolist()
                    default_wells = min(8, max(4, len(current_targets) // 3))
                else:
                    all_targets = st.session_state.expanded_df['target'].unique().tolist()
                    default_wells = min(8, max(4, len(all_targets) // 3))
                
                design_wells = st.number_input(
                    "设计孔数", 
                    min_value=1, 
                    max_value=50, 
                    value=default_wells, 
                    step=1,
                    help="组合方案必须严格按照该孔数分组"
                )
            
            # 开始合孔分析按钮
            if st.button("🚀 开始合孔分析", type="primary"):
                progress_container = st.container()
                with progress_container:
                    progress_bar = st.progress(0)
                    status_text = st.empty()
                    status_text.text("正在初始化计算...")
                try:
                    import time
                    start_time = time.time()
                    progress_bar.progress(10)
                    status_text.text("正在分析靶标组合...")
                    # 确保有勾选的序列
                    if not hasattr(st.session_state, 'sequences') or not st.session_state.sequences:
                        st.error("请先在上方选择要分析的序列")
                        progress_container.empty()
                        st.stop()
                    
                    # 过滤二聚体结果，只包含勾选序列的组合
                    selected_sequences = set(st.session_state.sequences)
                    filtered_dimer_results = [
                        result for result in st.session_state.sorted_results 
                        if result['seq1'] in selected_sequences and result['seq2'] in selected_sequences
                    ]
                    
                    st.session_state.optimized_combinations = optimize_combinations(
                        st.session_state.expanded_df, 
                        filtered_dimer_results,
                        max_combinations=10,
                        max_acceptable_deltaG=-3000.0,  # 固定值
                        min_acceptable_deltaG=min_acceptable_deltaG,
                        min_targets_per_well=min_targets_per_well,
                        max_targets_per_well=max_targets_per_well,
                        max_wells=design_wells
                    )
                    progress_bar.progress(100)
                    status_text.text("计算完成！")
                    end_time = time.time()
                    st.session_state.calculation_time = end_time - start_time
                    # calculation_count已经在optimize_combinations_with_ortools中设置
                    time.sleep(1)
                    progress_container.empty()
                except Exception as e:
                    st.error(f"计算过程中出现错误：{str(e)}")
                    progress_container.empty()
            
            # 显示结果
            if hasattr(st.session_state, 'optimized_combinations') and st.session_state.optimized_combinations:
                # 检查是否使用了智能分配算法
                if len(st.session_state.optimized_combinations) > 0 and 'total_wells' in st.session_state.optimized_combinations[0]:
                    # 智能分配算法结果
                    st.markdown("**📋 最优合孔方案**")
                    
                    # 显示计算统计信息
                    col_stats1, col_stats2 = st.columns(2)
                    with col_stats1:
                        if hasattr(st.session_state, 'calculation_time'):
                            st.info(f"⏱️ 计算耗时: {st.session_state.calculation_time:.2f} 秒")
                    with col_stats2:
                        if hasattr(st.session_state, 'calculation_count'):
                            st.info(f"🔢 计算次数: {st.session_state.calculation_count:,} 次")
                    
                    # 按分配方案分组
                    assignment_groups = {}
                    for combo in st.session_state.optimized_combinations:
                        assignment_id = combo.get('assignment_id', 0)
                        total_wells = combo.get('total_wells', 1)
                        assignment_valid = combo.get('assignment_valid', False)
                        
                        if assignment_id not in assignment_groups:
                            assignment_groups[assignment_id] = {
                                'wells': [],
                                'total_wells': total_wells,
                                'assignment_valid': assignment_valid,
                                'total_problematic': 0
                            }
                        
                        assignment_groups[assignment_id]['wells'].append(combo)
                        assignment_groups[assignment_id]['total_problematic'] += combo['problematic_count']
                    
                    # 创建方案概览表格
                    scheme_data = []
                    for assignment_id, assignment in assignment_groups.items():
                        status = "✅ 推荐" if assignment['assignment_valid'] else ("⚠️ 可考虑" if assignment['total_problematic'] <= 5 else "❌ 不推荐")
                        
                        # 计算总体问题deltaG
                        total_problematic_deltaG = 0
                        if assignment['wells']:
                            total_problematic_deltaG = assignment['wells'][0].get('total_problematic_deltaG', 0)
                        
                        scheme_data.append({
                            '方案': f"方案{assignment_id + 1}",
                            '孔数': assignment['total_wells'],
                            '状态': status,
                            '问题二聚体总数': assignment['total_problematic'],
                            '总体问题ΔG (cal/mol)': f"{total_problematic_deltaG:.1f}"
                        })
                    
                    if scheme_data:
                        st.dataframe(pd.DataFrame(scheme_data), use_container_width=True)
                        
                        # 显示总体deltaG说明
                        st.info("💡 **总体问题ΔG**: 所有孔内问题二聚体的ΔG总和，数值越大（越接近0）表示方案越优")
                    
                    # 显示最佳方案的详细信息
                    best_assignment = min(assignment_groups.items(), key=lambda x: x[1]['total_problematic'])
                    assignment_id, assignment = best_assignment
                    
                    st.markdown(f"**详细分配方案 (方案{assignment_id + 1}):**")
                    
                    # 显示该方案的总体deltaG
                    if assignment['wells']:
                        total_deltaG = assignment['wells'][0].get('total_problematic_deltaG', 0)
                        st.markdown(f"**🎯 该方案总体问题ΔG: {total_deltaG:.1f} cal/mol**")
                    
                    # 创建孔分配表格 - 去重显示
                    well_data = []
                    displayed_wells = set()
                    for well in assignment['wells']:
                        well_num = well.get('well_number', 1)
                        if well_num not in displayed_wells:
                            displayed_wells.add(well_num)
                            status = "✅" if well['problematic_count'] == 0 else ("⚠️" if well['problematic_count'] <= 2 else "❌")
                            
                            well_data.append({
                                '孔号': f"第{well_num}孔",
                                '靶标': ' + '.join([str(x) for x in well['target_names']]),
                                '序列数': len(well['sequences']),
                                '问题二聚体': well['problematic_count'],
                                '平均ΔG': f"{well['average']:.1f}",
                                '最差ΔG': f"{well['min']:.1f}",
                                '状态': status
                            })
                    
                    if well_data:
                        st.dataframe(pd.DataFrame(well_data), use_container_width=True)
                    
                    # 显示问题二聚体详情（如果有）
                    all_problematic = []
                    for well in assignment['wells']:
                        for dimer in well['problematic_dimers']:
                            all_problematic.append({
                                '孔号': f"第{well.get('well_number', 1)}孔",
                                '序列1': dimer['seq1'],
                                '序列2': dimer['seq2'],
                                'ΔG (cal/mol)': f"{dimer['deltaG']:.1f}",
                                '问题类型': dimer['problem_type']
                            })
                    
                    if all_problematic:
                        st.markdown("**⚠️ 问题二聚体详情:**")
                        st.dataframe(pd.DataFrame(all_problematic), use_container_width=True)
                
                else:
                    # 传统组合算法结果
                    st.markdown("**📋 最优组合方案**")
                    
                    # 创建组合概览表格
                    combo_data = []
                    for i, combo in enumerate(st.session_state.optimized_combinations[:max_combinations_to_show], 1):
                        status = "✅ 推荐" if combo['problematic_count'] == 0 else ("⚠️ 可考虑" if combo['problematic_count'] <= 2 else "❌ 不推荐")
                        combo_data.append({
                            '方案': f"组合{i}",
                            '靶标': ' + '.join([str(x) for x in combo['target_names']]),
                            '序列数': len(combo['sequences']),
                            '问题二聚体': combo['problematic_count'],
                            '平均ΔG': f"{combo['average']:.1f}",
                            '最差ΔG': f"{combo['min']:.1f}",
                            '状态': status
                        })
                    
                    if combo_data:
                        st.dataframe(pd.DataFrame(combo_data), use_container_width=True)
                    
                    # 显示最佳组合的问题二聚体详情
                    if st.session_state.optimized_combinations and st.session_state.optimized_combinations[0]['problematic_count'] > 0:
                        best_combo = st.session_state.optimized_combinations[0]
                        st.markdown("**⚠️ 问题二聚体详情:**")
                        
                        problem_data = []
                        for dimer in best_combo['problematic_dimers']:
                            problem_data.append({
                                '序列1': dimer['seq1'],
                                '序列2': dimer['seq2'],
                                'ΔG (cal/mol)': f"{dimer['deltaG']:.1f}",
                                '问题类型': dimer['problem_type']
                            })
                        
                        if problem_data:
                            st.dataframe(pd.DataFrame(problem_data), use_container_width=True)
                
                # 添加下载组合建议的功能
                if len(st.session_state.optimized_combinations) > 0:
                    # 创建用于下载的DataFrame
                    combo_download_data = []
                    for i, combo in enumerate(st.session_state.optimized_combinations, 1):
                        combo_download_data.append({
                            '组合编号': i,
                            '靶标组合': ' + '.join([str(x) for x in combo['target_names']]),
                            '包含序列': ', '.join(combo['sequences']),
                            '二聚体数量': combo['count'],
                            '总ΔG (cal/mol)': round(combo['total'], 2),
                            '平均ΔG (cal/mol)': round(combo['average'], 2),
                            '最小ΔG (cal/mol)': round(combo['min'], 2),
                            '最大ΔG (cal/mol)': round(combo['max'], 2)
                        })
                    
                    combo_download_df = pd.DataFrame(combo_download_data)
                    
                    # 创建Excel文件的缓冲区
                    buffer = io.BytesIO()
                    with pd.ExcelWriter(buffer, engine='openpyxl') as writer:
                        combo_download_df.to_excel(writer, index=False, sheet_name='组合优化建议')
                        
                        # 获取工作表对象以调整列宽
                        worksheet = writer.sheets['组合优化建议']
                        
                        # 自动调整列宽
                        for column in worksheet.columns:
                            max_length = 0
                            column_letter = column[0].column_letter
                            for cell in column:
                                try:
                                    if len(str(cell.value)) > max_length:
                                        max_length = len(str(cell.value))
                                except:
                                    pass
                            adjusted_width = min(max_length + 2, 50)
                            worksheet.column_dimensions[column_letter].width = adjusted_width
                    
                    buffer.seek(0)
                    
                    # 生成下载文件名
                    from datetime import datetime
                    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                    combo_download_filename = f"靶标组合优化建议_{timestamp}.xlsx"
                    
                    st.download_button(
                        label="📥 下载组合优化建议 (Excel格式)",
                        data=buffer.getvalue(),
                        file_name=combo_download_filename,
                        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                        help="下载所有组合优化建议"
                    )
            else:
                st.info("未找到合适的组合建议，可能是因为靶标数量不足或没有检测到二聚体。")
            
            st.markdown("---")
        
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
        
        # 添加下载分析结果的功能
        if len(filtered_results) > 0:
            col_download1, col_download2 = st.columns(2)
            
            with col_download1:
                # 创建用于下载的DataFrame
                download_data = []
                for res in filtered_results:
                    download_data.append({
                        '序列1': res['seq1'],
                        '序列2': res['seq2'],
                        'Tm (°C)': round(res['Tm'], 2),
                        'ΔG (cal/mol)': round(res['deltaG'], 2),
                        '二聚体结构': res['结构']
                    })
                
                download_df = pd.DataFrame(download_data)
                
                # 创建Excel文件的缓冲区
                buffer = io.BytesIO()
                with pd.ExcelWriter(buffer, engine='openpyxl') as writer:
                    download_df.to_excel(writer, index=False, sheet_name='二聚体分析结果')
                    
                    # 获取工作表对象以调整列宽
                    worksheet = writer.sheets['二聚体分析结果']
                    
                    # 自动调整列宽
                    for column in worksheet.columns:
                        max_length = 0
                        column_letter = column[0].column_letter
                        for cell in column:
                            try:
                                if len(str(cell.value)) > max_length:
                                    max_length = len(str(cell.value))
                            except:
                                pass
                        adjusted_width = min(max_length + 2, 50)  # 限制最大宽度为50
                        worksheet.column_dimensions[column_letter].width = adjusted_width
                
                buffer.seek(0)
                
                # 生成下载文件名
                from datetime import datetime
                timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                download_filename = f"二聚体分析结果_{timestamp}.xlsx"
                
                st.download_button(
                    label="📥 下载二聚体分析结果 (Excel)",
                    data=buffer.getvalue(),
                    file_name=download_filename,
                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                    help="下载当前筛选条件下的所有二聚体分析结果"
                )
            
            with col_download2:
                # 创建CSV格式的下载
                csv_data = download_df.to_csv(index=False, encoding='utf-8-sig')
                csv_filename = f"二聚体分析结果_{timestamp}.csv"
                
                st.download_button(
                    label="📥 下载二聚体分析结果 (CSV)",
                    data=csv_data.encode('utf-8-sig'),
                    file_name=csv_filename,
                    mime="text/csv",
                    help="下载当前筛选条件下的所有二聚体分析结果（CSV格式）"
                )
        
        # 显示选项控制
        col_display1, col_display2 = st.columns([2, 1])
        with col_display1:
            display_count = st.selectbox(
                "显示结果数量:",
                [50, 100, 200, 500, "全部"],
                index=0,
                help="选择要显示的二聚体结果数量，显示过多可能导致页面卡顿"
            )
        
        with col_display2:
            if len(filtered_results) > 50:
                st.info(f"💡 为避免卡顿，默认显示前50个结果")
        
        # 确定要显示的结果数量
        if display_count == "全部":
            results_to_show = filtered_results
        else:
            results_to_show = filtered_results[:display_count]
        
        # 显示当前显示的结果数量
        if len(results_to_show) < len(filtered_results):
            st.write(f"当前显示: {len(results_to_show)} / {len(filtered_results)} 个结果")
        
        # 显示筛选后的结果
        for res in results_to_show:
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

st.markdown("---")
