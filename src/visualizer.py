#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
可视化模块
Visualization Module

功能: 生成工作流结果的可视化图表
"""

import os
import sys
import pandas as pd
import numpy as np

plt = None
sns = None

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import seaborn as sns
except ImportError:
    pass


def generate_plots(final_df, qsar_results=None, output_path='results/evaluation_summary.png'):
    """
    生成可视化图表
    
    参数:
        final_df: 最终结果数据框
        qsar_results: QSAR结果数据框
        output_path: 输出图片路径
    """
    print("=" * 60)
    print("步骤6: 结果可视化")
    print("=" * 60)
    
    if final_df is None or len(final_df) == 0:
        print("警告: 没有数据可供可视化")
        return
    
    if plt is None:
        print("警告: 缺少matplotlib库，跳过可视化步骤")
        print("请安装: pip install matplotlib seaborn")
        return
    
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    try:
        _plot_ic50_distribution(final_df, axes[0, 0])
    except Exception as e:
        print(f"IC50分布图生成失败: {e}")
        axes[0, 0].text(0.5, 0.5, 'IC50 Distribution\n(No Data)', ha='center', va='center')
    
    try:
        _plot_chemical_space(final_df, axes[0, 1])
    except Exception as e:
        print(f"化学空间图生成失败: {e}")
        axes[0, 1].text(0.5, 0.5, 'Chemical Space\n(No Data)', ha='center', va='center')
    
    try:
        _plot_sa_score_distribution(final_df, axes[1, 0])
    except Exception as e:
        print(f"SA Score分布图生成失败: {e}")
        axes[1, 0].text(0.5, 0.5, 'SA Score Distribution\n(No Data)', ha='center', va='center')
    
    try:
        _plot_top_candidates(final_df, axes[1, 1])
    except Exception as e:
        print(f"候选分子图生成失败: {e}")
        axes[1, 1].text(0.5, 0.5, 'Top Candidates\n(No Data)', ha='center', va='center')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"可视化图表已保存至: {output_path}")


def _plot_ic50_distribution(df, ax):
    """绘制IC50分布图"""
    if 'pred_IC50' not in df.columns:
        ax.text(0.5, 0.5, 'No IC50 Data', ha='center', va='center')
        return
    
    sns.histplot(data=df, x='pred_IC50', bins=20, color='skyblue', kde=True, ax=ax)
    ax.set_xlabel('Predicted IC50 (nM)')
    ax.set_ylabel('Count')
    ax.set_title('IC50 Distribution')
    ax.axvline(x=100, color='red', linestyle='--', label='Threshold (100 nM)')
    ax.legend()


def _plot_chemical_space(df, ax):
    """绘制化学空间分布图"""
    if 'MW' not in df.columns or 'LogP' not in df.columns:
        if 'SMILES' in df.columns:
            df = _calculate_properties(df)
        else:
            ax.text(0.5, 0.5, 'No Property Data', ha='center', va='center')
            return
    
    hue_col = 'sa_score' if 'sa_score' in df.columns else None
    sns.scatterplot(data=df, x='MW', y='LogP', hue=hue_col, 
                    palette='viridis', alpha=0.6, ax=ax)
    ax.set_xlabel('Molecular Weight (Da)')
    ax.set_ylabel('LogP')
    ax.set_title('Chemical Space Distribution')
    if hue_col:
        ax.legend(title='SA Score')


def _plot_sa_score_distribution(df, ax):
    """绘制SA Score分布图"""
    if 'sa_score' not in df.columns:
        ax.text(0.5, 0.5, 'No SA Score Data', ha='center', va='center')
        return
    
    sns.boxplot(data=df, y='sa_score', color='lightgreen', ax=ax)
    ax.set_ylabel('SA Score')
    ax.set_title('SA Score Distribution')
    ax.axhline(y=6.0, color='red', linestyle='--', label='Threshold (6.0)')
    ax.legend()


def _plot_top_candidates(df, ax):
    """绘制Top候选分子"""
    if 'pred_IC50' not in df.columns:
        ax.text(0.5, 0.5, 'No IC50 Data', ha='center', va='center')
        return
    
    top_n = min(20, len(df))
    if 'sa_score' in df.columns:
        top_df = df.nsmallest(top_n, 'sa_score')
    else:
        top_df = df.nsmallest(top_n, 'pred_IC50')
    
    top_df = top_df.reset_index(drop=True)
    top_df['rank'] = range(1, len(top_df) + 1)
    
    sns.barplot(data=top_df, x='rank', y='pred_IC50', color='steelblue', ax=ax)
    ax.set_xlabel('Rank')
    ax.set_ylabel('Predicted IC50 (nM)')
    ax.set_title(f'Top {top_n} Candidates')
    ax.tick_params(axis='x', rotation=45)


def _calculate_properties(df):
    """计算分子性质"""
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors
        
        mw_list = []
        logp_list = []
        
        for smiles in df['SMILES']:
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    mw_list.append(Descriptors.MolWt(mol))
                    logp_list.append(Descriptors.MolLogP(mol))
                else:
                    mw_list.append(None)
                    logp_list.append(None)
            except:
                mw_list.append(None)
                logp_list.append(None)
        
        df['MW'] = mw_list
        df['LogP'] = logp_list
        
    except ImportError:
        pass
    
    return df


def generate_summary_report(df, output_path='results/summary_report.txt'):
    """
    生成文本摘要报告
    
    参数:
        df: 结果数据框
        output_path: 输出文件路径
    """
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write("=" * 60 + "\n")
        f.write("分子筛选工作流结果摘要\n")
        f.write("=" * 60 + "\n\n")
        
        f.write(f"总分子数: {len(df)}\n\n")
        
        if 'pred_IC50' in df.columns:
            f.write("IC50统计:\n")
            f.write(f"  最小值: {df['pred_IC50'].min():.2f} nM\n")
            f.write(f"  最大值: {df['pred_IC50'].max():.2f} nM\n")
            f.write(f"  平均值: {df['pred_IC50'].mean():.2f} nM\n")
            f.write(f"  中位数: {df['pred_IC50'].median():.2f} nM\n\n")
        
        if 'sa_score' in df.columns:
            f.write("SA Score统计:\n")
            f.write(f"  最小值: {df['sa_score'].min():.2f}\n")
            f.write(f"  最大值: {df['sa_score'].max():.2f}\n")
            f.write(f"  平均值: {df['sa_score'].mean():.2f}\n\n")
        
        if 'binding_affinity' in df.columns:
            f.write("结合能统计:\n")
            f.write(f"  最佳: {df['binding_affinity'].min():.2f} kcal/mol\n")
            f.write(f"  平均: {df['binding_affinity'].mean():.2f} kcal/mol\n\n")
        
        f.write("=" * 60 + "\n")
    
    print(f"摘要报告已保存至: {output_path}")


if __name__ == "__main__":
    test_df = pd.DataFrame({
        'SMILES': ['CCO', 'CCCO', 'CCCCO'],
        'pred_IC50': [45.2, 82.1, 35.8],
        'sa_score': [2.1, 3.5, 2.8]
    })
    generate_plots(test_df, output_path='test_visualization.png')
    print("测试可视化完成")
