import pandas as pd
import numpy as np
import os
import sys
import logging
import argparse
from typing import List, Tuple
from datetime import datetime

# 创建全局logger对象
logger = None

def setup_logging(importance_thresholds: str, occurrence_thresholds: str) -> logging.Logger:
    """
    设置日志配置
    
    Args:
        importance_thresholds (str): 重要性阈值字符串
        occurrence_thresholds (str): 出现频率阈值字符串
    
    Returns:
        logging.Logger: 配置好的logger对象
    """
    # 创建日志目录
    os.makedirs('logs', exist_ok=True)
    
    # 生成日志文件名，包含时间戳和参数信息
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_filename = f'logs/merge_network_{timestamp}.log'
    
    # 配置日志格式
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s: %(message)s',
        handlers=[
            logging.FileHandler(log_filename, mode='w', encoding='utf-8'),
            logging.StreamHandler(sys.stdout)
        ]
    )
    
    # 获取logger
    global logger
    logger = logging.getLogger(__name__)
    
    # 记录脚本启动信息
    logger.info("=" * 80)
    logger.info("网络合并分析启动")
    logger.info("=" * 80)
    logger.info("参数配置:")
    logger.info(f"  重要性阈值: {importance_thresholds}")
    logger.info(f"  出现频率阈值: {occurrence_thresholds}")
    logger.info("-" * 80)
    
    return logger

def parse_args():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(description='合并GRN网络文件')
    parser.add_argument('--num-runs', type=int, default=20,
                      help='GRN运行次数')
    parser.add_argument('--importance-thresholds', type=str, 
                      default="0.1,0.3,0.5,0.7,1.0",
                      help='重要性阈值，用逗号分隔')
    parser.add_argument('--occurrence-thresholds', type=str, 
                      default="5,10,15,20",
                      help='出现频率阈值，用逗号分隔')
    return parser.parse_args()

def smart_read_network_file(file_path: str) -> pd.DataFrame:
    """
    智能读取网络文件，保持原始格式
    
    Args:
        file_path (str): 网络文件路径
    
    Returns:
        pd.DataFrame: 处理后的网络DataFrame
    """
    try:
        # 读取文件，保持原始列名
        df = pd.read_csv(
            file_path, 
            sep='\t', 
            dtype={
                'TF': str, 
                'target': str, 
                'importance': np.float32
            }
        )
        
        logger.info(f"成功读取文件: {file_path}")
        logger.info(f"  文件预览: \n{df.head()}")
        logger.info(f"  列名: {df.columns.tolist()}")
        
        return df
    
    except Exception as e:
        logger.error(f"读取文件 {file_path} 时发生错误: {e}")
        raise

def load_network_files(num_runs: int) -> List[pd.DataFrame]:
    """
    加载并验证指定数量的网络文件
    
    Args:
        num_runs (int): 运行次数
    
    Returns:
        List[pd.DataFrame]: 有效的网络DataFrame列表
    """
    df_list = []
    valid_files = 0
    
    for run in range(1, num_runs + 1):
        file_path = f"GRN_output/raw_grn/run_{run}_s1_adj.tsv"
        
        try:
            # 使用智能读取函数
            df = smart_read_network_file(file_path)
            
            # 数据验证和过滤
            valid_df = df[
                (df['importance'] > 0) &  # 大于0的重要性值
                (df['TF'].notna()) & 
                (df['target'].notna()) & 
                (df['TF'] != df['target'])  # 排除自环
            ]
            
            if len(valid_df) > 0:
                df_list.append(valid_df)
                valid_files += 1
                
                logger.info(f"成功加载文件: {file_path}")
                logger.info(f"  总行数: {len(df)}")
                logger.info(f"  有效行数: {len(valid_df)}")
                logger.info(f"  重要性范围: {valid_df['importance'].min():.4f} - {valid_df['importance'].max():.4f}")
            else:
                logger.warning(f"文件 {file_path} 无有效数据")
        
        except Exception as e:
            logger.error(f"加载文件 {file_path} 时发生错误: {e}")
    
    # 检查有效文件数量
    if valid_files < num_runs * 0.8:
        raise ValueError(f"仅加载了 {valid_files}/{num_runs} 个有效文件，低于预期80%")
    
    return df_list

def merge_networks(
    df_list: List[pd.DataFrame], 
    importance_thresholds: List[float], 
    occurrence_thresholds: List[int]
):
    """
    合并网络并根据重要性和出现频率生成多个网络文件
    
    Args:
        df_list (List[pd.DataFrame]): 网络DataFrame列表
        importance_thresholds (List[float]): 重要性阈值列表
        occurrence_thresholds (List[int]): 出现频率阈值列表
    """
    # 合并所有网络
    full_network = pd.concat(df_list)
    
    # 计算初始统计信息
    initial_tfs = full_network['TF'].nunique()
    initial_targets = full_network['target'].nunique()
    initial_edges = len(full_network)
    
    logger.info("初始网络统计:")
    logger.info(f"  转录因子数量: {initial_tfs}")
    logger.info(f"  靶基因数量: {initial_targets}")
    logger.info(f"  调控关系数量: {initial_edges}")
    
    # 计算网络统计
    network_stats = full_network.groupby(['TF', 'target']).agg(
        Mean_Importance=('importance', 'mean'),
        Std_Importance=('importance', 'std'),
        Occurrence=('importance', 'size')
    ).reset_index()
    
    # 保存全局网络统计
    os.makedirs('GRN_output/stats_report', exist_ok=True)
    network_stats.to_csv('GRN_output/stats_report/full_network_stats.tsv', sep='\t', index=False)
    
    # 创建统计报告文件
    stats_report_path = 'GRN_output/stats_report/network_filtering_stats.tsv'
    with open(stats_report_path, 'w') as f:
        f.write("Importance_Threshold\tOccurrence_Threshold\tTF_Count\tTarget_Count\tEdge_Count\tMean_Importance\tMin_Importance\tMax_Importance\n")
    
    logger.info(f"全局网络统计:")
    logger.info(f"  总网络边数: {len(network_stats)}")
    logger.info(f"  平均重要性: {network_stats['Mean_Importance'].mean():.4f}")
    logger.info(f"  重要性范围: {network_stats['Mean_Importance'].min():.4f} - {network_stats['Mean_Importance'].max():.4f}")
    
    # 根据阈值生成网络文件
    for imp_threshold in importance_thresholds:
        for occ_threshold in occurrence_thresholds:
            logger.info(f"\n处理阈值组合: imp={imp_threshold:.3f}, occ={occ_threshold}")
            
            # 过滤网络
            filtered_network = network_stats[
                (network_stats['Mean_Importance'] >= imp_threshold) & 
                (network_stats['Occurrence'] >= occ_threshold)
            ]
            
            # 创建输出目录
            output_dir = f"GRN_output/final_results/imp_{imp_threshold:.3f}_occ_{occ_threshold}"
            os.makedirs(output_dir, exist_ok=True)
            
            # 输出网络文件
            if len(filtered_network) > 0:
                output_path = f"{output_dir}/network.tsv"
                
                # 保留原始列名和格式
                output_network = filtered_network[['TF', 'target', 'Mean_Importance']].copy()
                output_network.columns = ['TF', 'target', 'importance']
                output_network['importance'] = output_network['importance'].round(4)
                
                # 计算统计信息
                tf_count = output_network['TF'].nunique()
                target_count = output_network['target'].nunique()
                edge_count = len(output_network)
                mean_imp = output_network['importance'].mean()
                min_imp = output_network['importance'].min()
                max_imp = output_network['importance'].max()
                
                # 记录详细统计信息
                logger.info(f"  网络统计:")
                logger.info(f"    转录因子数量: {tf_count} ({(tf_count/initial_tfs*100):.1f}% of initial)")
                logger.info(f"    靶基因数量: {target_count} ({(target_count/initial_targets*100):.1f}% of initial)")
                logger.info(f"    调控关系数量: {edge_count} ({(edge_count/initial_edges*100):.1f}% of initial)")
                logger.info(f"    重要性分布:")
                logger.info(f"      平均值: {mean_imp:.4f}")
                logger.info(f"      最小值: {min_imp:.4f}")
                logger.info(f"      最大值: {max_imp:.4f}")
                
                # 保存统计信息到报告文件
                with open(stats_report_path, 'a') as f:
                    f.write(f"{imp_threshold}\t{occ_threshold}\t{tf_count}\t{target_count}\t{edge_count}\t{mean_imp:.4f}\t{min_imp:.4f}\t{max_imp:.4f}\n")
                
                # 保存网络文件
                output_network.to_csv(
                    output_path, 
                    sep='\t', 
                    index=False
                )
                
                # 保存该阈值组合下的TF和target列表
                tf_list_path = f"{output_dir}/transcription_factors.txt"
                target_list_path = f"{output_dir}/target_genes.txt"
                
                with open(tf_list_path, 'w') as f:
                    f.write('\n'.join(sorted(output_network['TF'].unique())))
                
                with open(target_list_path, 'w') as f:
                    f.write('\n'.join(sorted(output_network['target'].unique())))
                
                logger.info(f"  输出文件:")
                logger.info(f"    网络文件: {output_path}")
                logger.info(f"    TF列表: {tf_list_path}")
                logger.info(f"    靶基因列表: {target_list_path}")
            else:
                logger.warning(f"  阈值 imp={imp_threshold}, occ={occ_threshold} 未生成网络")
    
    logger.info("\n网络过滤统计报告已保存至: " + stats_report_path)

def parse_thresholds(threshold_str: str) -> List[float]:
    """
    解析阈值字符串，处理空格并转换为浮点数列表
    
    Args:
        threshold_str (str): 逗号分隔的阈值字符串
    
    Returns:
        List[float]: 浮点数阈值列表
    """
    try:
        # 分割字符串，去除空格，过滤空字符串
        values = [x.strip() for x in threshold_str.split(',') if x.strip()]
        # 转换为浮点数
        return [float(x) for x in values]
    except ValueError as e:
        logger.error(f"阈值解析错误，输入字符串: '{threshold_str}'")
        raise ValueError(f"无法解析阈值字符串: {e}")

def main():
    # 解析命令行参数
    args = parse_args()
    
    try:
        # 设置日志
        setup_logging(args.importance_thresholds, args.occurrence_thresholds)
        
        # 解析阈值列表
        importance_thresholds = parse_thresholds(args.importance_thresholds)
        occurrence_thresholds = [int(x.strip()) for x in args.occurrence_thresholds.split(',') if x.strip()]
        
        logger.info(f"开始网络合并分析:")
        logger.info(f"  运行次数: {args.num_runs}")
        logger.info(f"  重要性阈值: {importance_thresholds}")
        logger.info(f"  出现频率阈值: {occurrence_thresholds}")
        
        # 记录系统信息
        logger.info("\n系统信息:")
        logger.info(f"  Python版本: {sys.version}")
        logger.info(f"  工作目录: {os.getcwd()}")
        logger.info(f"  命令行参数: {' '.join(sys.argv)}")
        
        # 加载网络文件
        logger.info("\n开始加载网络文件...")
        network_dfs = load_network_files(args.num_runs)
        logger.info(f"成功加载 {len(network_dfs)} 个网络文件")
        
        # 合并网络
        logger.info("\n开始合并网络...")
        merge_networks(
            network_dfs, 
            importance_thresholds, 
            occurrence_thresholds
        )
        
        logger.info("\n网络合并分析完成！")
        logger.info("=" * 80)
    
    except Exception as e:
        logger.error(f"网络合并分析失败: {e}")
        logger.exception("详细错误信息:")
        sys.exit(1)

if __name__ == '__main__':
    main()
