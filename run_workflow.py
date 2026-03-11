#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
分子筛选工作流统一入口程序
Molecular Screening Workflow - Unified Entry Point

作者: Molecular Screening Team
版本: 1.0.0
日期: 2026-03-11

功能说明:
    本脚本整合了分子筛选工作流的所有步骤，包括：
    1. RNN分子生成
    2. QSAR活性预测
    3. ADMET性质过滤
    4. 分子对接
    5. SA Score评估
    6. 结果可视化

使用方法:
    python run_workflow.py --config config/workflow_config.yaml
    python run_workflow.py --demo  # 运行演示模式
"""

import os
import sys
import time
import logging
import argparse
import traceback
from datetime import datetime
from pathlib import Path

# 添加src目录到系统路径
src_dir = os.path.join(os.path.dirname(__file__), 'src')
sys.path.insert(0, src_dir)

# 导入工作流模块
try:
    from rnn_workflow import run_rnn_generation
    from qsar_engine import qsar_predict_and_filter
    from admet_engine import run_admet_filter
    from submit_job_docking import submit_docking_job
    from sa_scorer import run_sa_filter
    from visualizer import generate_plots
except ImportError as e:
    print(f"错误: 无法导入工作流模块 - {e}")
    print("请确保所有依赖模块已正确安装")
    sys.exit(1)


class WorkflowLogger:
    """工作流日志记录器"""
    
    def __init__(self, log_dir='logs'):
        self.log_dir = log_dir
        os.makedirs(log_dir, exist_ok=True)
        
        # 创建日志文件名
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        self.log_file = os.path.join(log_dir, f'workflow_{timestamp}.log')
        
        # 配置日志
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(self.log_file, encoding='utf-8'),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger(__name__)
        
    def log_step(self, step_name, status, duration=None, output_path=None):
        """记录工作流步骤"""
        msg = f"[{step_name}] - {status}"
        if duration:
            msg += f" (耗时: {duration:.2f}秒)"
        if output_path:
            msg += f" -> {output_path}"
        self.logger.info(msg)


class MolecularScreeningWorkflow:
    """分子筛选工作流主类"""
    
    def __init__(self, config, demo_mode=False):
        self.config = config
        self.demo_mode = demo_mode
        self.logger = WorkflowLogger()
        self.results = {}
        
        # 初始化目录
        self._init_directories()
        
    def _init_directories(self):
        """初始化工作流目录"""
        dirs = [
            'results',
            'results/QSAR',
            'results/ADMET', 
            'results/docking',
            'results/SA',
            'logs',
            'models'
        ]
        for d in dirs:
            os.makedirs(d, exist_ok=True)
        self.logger.log_step("初始化", "目录结构创建完成")
    
    def run_step(self, step_func, step_name, *args, **kwargs):
        """运行单个工作流步骤"""
        self.logger.log_step(step_name, "开始执行")
        start_time = time.time()
        
        try:
            result = step_func(*args, **kwargs)
            duration = time.time() - start_time
            self.logger.log_step(step_name, "执行成功", duration)
            return result
        except Exception as e:
            duration = time.time() - start_time
            error_msg = f"执行失败: {str(e)}"
            self.logger.log_step(step_name, error_msg, duration)
            self.logger.logger.error(traceback.format_exc())
            raise
    
    def run_rnn_generation(self):
        """步骤1: RNN分子生成"""
        self.logger.log_step("RNN分子生成", "开始")
        
        start_time = time.time()
        
        try:
            result = run_rnn_generation(self.config, demo_mode=self.demo_mode)
            duration = time.time() - start_time
            output_path = 'results/generated_smiles_from_1994mols.csv'
            self.logger.log_step("RNN分子生成", "完成", duration, output_path)
            self.results['rnn_output'] = output_path
            return result
        except Exception as e:
            duration = time.time() - start_time
            self.logger.log_step("RNN分子生成", f"失败: {str(e)}", duration)
            raise
    
    def run_qsar_prediction(self):
        """步骤2: QSAR活性预测"""
        self.logger.log_step("QSAR活性预测", "开始")
        
        start_time = time.time()
        
        try:
            train_data = self.config.get('qsar_train_data', 'data/ChEMBL_Cleaned_akt1-activities_zip.csv')
            predict_data = self.config.get('qsar_predict_data', 'results/generated_smiles_from_1994mols.csv')
            model_path = self.config.get('qsar_model_path', 'models/ultimate_ensemble_qsar_model.pkl')
            
            result_df = qsar_predict_and_filter(train_data, predict_data, model_path, demo_mode=self.demo_mode)
            
            duration = time.time() - start_time
            output_path = 'results/QSAR/top_100_generated_molecules_with_predictions.csv'
            self.logger.log_step("QSAR活性预测", "完成", duration, output_path)
            self.results['qsar_output'] = output_path
            return result_df
        except Exception as e:
            duration = time.time() - start_time
            self.logger.log_step("QSAR活性预测", f"失败: {str(e)}", duration)
            raise
    
    def run_admet_filter(self):
        """步骤3: ADMET性质过滤"""
        self.logger.log_step("ADMET性质过滤", "开始")
        
        start_time = time.time()
        
        try:
            input_df = self.results.get('qsar_df')
            if input_df is None:
                input_path = 'results/QSAR/top_100_generated_molecules_with_predictions.csv'
                import pandas as pd
                input_df = pd.read_csv(input_path)
            
            admet_df = run_admet_filter(input_df)
            
            duration = time.time() - start_time
            output_path = 'results/ADMET/top_100_generated_molecules_only_admet_passed.csv'
            self.logger.log_step("ADMET性质过滤", "完成", duration, output_path)
            self.results['admet_output'] = output_path
            return admet_df
        except Exception as e:
            duration = time.time() - start_time
            self.logger.log_step("ADMET性质过滤", f"失败: {str(e)}", duration)
            raise
    
    def run_docking(self):
        """步骤4: 分子对接"""
        self.logger.log_step("分子对接", "开始")
        
        start_time = time.time()
        
        try:
            docking_engine = self.config.get('docking_engine', ['vina'])
            result = submit_docking_job(docking_engine, demo_mode=self.demo_mode)
            
            duration = time.time() - start_time
            self.logger.log_step("分子对接", "任务完成", duration)
            self.results['docking_result'] = result
            return result
        except Exception as e:
            duration = time.time() - start_time
            self.logger.log_step("分子对接", f"失败: {str(e)}", duration)
            self.logger.logger.warning("分子对接步骤失败，继续执行后续步骤")
            return None
    
    def run_sa_scoring(self):
        """步骤5: SA Score评估"""
        self.logger.log_step("SA Score评估", "开始")
        
        start_time = time.time()
        
        try:
            import pandas as pd
            input_path = 'results/docking/docking_result.csv'
            
            if not os.path.exists(input_path):
                self.logger.logger.warning("对接结果文件不存在，使用ADMET结果")
                input_path = self.results.get('admet_output', 'results/ADMET/top_100_generated_molecules_only_admet_passed.csv')
            
            df = pd.read_csv(input_path)
            sa_df = run_sa_filter(df)
            
            duration = time.time() - start_time
            output_path = 'results/SA/sa_filtered.csv'
            self.logger.log_step("SA Score评估", "完成", duration, output_path)
            self.results['sa_output'] = output_path
            return sa_df
        except Exception as e:
            duration = time.time() - start_time
            self.logger.log_step("SA Score评估", f"失败: {str(e)}", duration)
            raise
    
    def run_visualization(self):
        """步骤6: 结果可视化"""
        self.logger.log_step("结果可视化", "开始")
        
        start_time = time.time()
        
        try:
            final_df = self.results.get('sa_df')
            qsar_results = self.results.get('qsar_df')
            
            generate_plots(final_df, qsar_results)
            
            duration = time.time() - start_time
            output_path = 'results/evaluation_summary.png'
            self.logger.log_step("结果可视化", "完成", duration, output_path)
            return output_path
        except Exception as e:
            duration = time.time() - start_time
            self.logger.log_step("结果可视化", f"失败: {str(e)}", duration)
            self.logger.logger.warning("可视化步骤失败，但工作流已完成")
    
    def run_full_workflow(self):
        """运行完整工作流"""
        self.logger.logger.info("=" * 60)
        self.logger.logger.info("分子筛选工作流开始执行")
        self.logger.logger.info("=" * 60)
        
        workflow_start = time.time()
        
        try:
            # 步骤1: RNN分子生成
            self.run_rnn_generation()
            
            # 步骤2: QSAR活性预测
            qsar_df = self.run_qsar_prediction()
            self.results['qsar_df'] = qsar_df
            
            # 步骤3: ADMET性质过滤
            admet_df = self.run_admet_filter()
            self.results['admet_df'] = admet_df
            
            # 步骤4: 分子对接
            self.run_docking()
            
            # 步骤5: SA Score评估
            sa_df = self.run_sa_scoring()
            self.results['sa_df'] = sa_df
            
            # 步骤6: 结果可视化
            self.run_visualization()
            
            workflow_duration = time.time() - workflow_start
            
            self.logger.logger.info("=" * 60)
            self.logger.logger.info(f"工作流执行完成！总耗时: {workflow_duration:.2f}秒")
            self.logger.logger.info("=" * 60)
            
            # 输出结果摘要
            self._print_summary()
            
            return self.results
            
        except Exception as e:
            workflow_duration = time.time() - workflow_start
            self.logger.logger.error("=" * 60)
            self.logger.logger.error(f"工作流执行失败: {str(e)}")
            self.logger.logger.error(f"总耗时: {workflow_duration:.2f}秒")
            self.logger.logger.error("=" * 60)
            raise
    
    def _print_summary(self):
        """打印工作流结果摘要"""
        self.logger.logger.info("\n结果摘要:")
        self.logger.logger.info("-" * 40)
        
        for key, value in self.results.items():
            if isinstance(value, str) and os.path.exists(value):
                self.logger.logger.info(f"{key}: {value}")
        
        self.logger.logger.info(f"\n日志文件: {self.logger.log_file}")
        self.logger.logger.info("工作流执行完成！")


def load_config(config_path):
    """加载配置文件"""
    try:
        import yaml
        has_yaml = True
    except ImportError:
        has_yaml = False
    
    if not os.path.exists(config_path):
        print(f"警告: 配置文件 {config_path} 不存在，使用默认配置")
        return {}
    
    if not has_yaml:
        print("警告: 缺少pyyaml模块，使用默认配置")
        print("请安装: pip install pyyaml")
        return {}
    
    with open(config_path, 'r', encoding='utf-8') as f:
        config = yaml.safe_load(f)
    
    return config or {}


def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description='分子筛选工作流 - 统一入口程序',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例用法:
  # 使用默认配置运行完整工作流
  python run_workflow.py
  
  # 使用自定义配置文件
  python run_workflow.py --config config/workflow_config.yaml
  
  # 运行演示模式
  python run_workflow.py --demo
  
  # 只运行特定步骤
  python run_workflow.py --steps rnn qsar admet
        """
    )
    
    parser.add_argument(
        '--config',
        type=str,
        default='config/workflow_config.yaml',
        help='配置文件路径 (默认: config/workflow_config.yaml)'
    )
    
    parser.add_argument(
        '--demo',
        action='store_true',
        help='运行演示模式（使用最小数据集）'
    )
    
    parser.add_argument(
        '--steps',
        nargs='+',
        choices=['rnn', 'qsar', 'admet', 'docking', 'sa', 'viz'],
        help='只运行指定步骤'
    )
    
    args = parser.parse_args()
    
    # 加载配置
    config = load_config(args.config)
    
    # 演示模式配置
    demo_mode = args.demo
    if demo_mode:
        config['demo_mode'] = True
        config['qsar_train_data'] = 'demo/data/sample_training.csv'
        config['qsar_predict_data'] = 'demo/data/sample_molecules.csv'
        print("运行演示模式...")
    
    # 创建工作流实例
    workflow = MolecularScreeningWorkflow(config, demo_mode=demo_mode)
    
    # 运行工作流
    try:
        results = workflow.run_full_workflow()
        print("\n✅ 工作流执行成功！")
        print(f"查看详细日志: {workflow.logger.log_file}")
        return 0
    except Exception as e:
        print(f"\n❌ 工作流执行失败: {str(e)}")
        print(f"查看错误日志: {workflow.logger.log_file}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
