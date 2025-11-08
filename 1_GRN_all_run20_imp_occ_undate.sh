#!/bin/bash
### —— 用户参数配置 ###
export NUM_RUNS=20
# export THREADS=10  # 删除这行，使用新的线程配置

# 修改阈值配置
export IMPORTANCE_THRESHOLDS="0.1,0.3,0.5,0.7,1.0"
export OCCURRENCE_THRESHOLDS="5,10,15,20"

### —— 系统资源配置 —— ###
# 基础配置
export PYTHONHASHSEED=42
export HDF5_USE_FILE_LOCKING=FALSE

# 内存限制配置（限制在150GB以内）
export MAX_MEMORY_USAGE=150          # 最大内存使用限制（GB）
export JAVA_OPTS="-Xmx40g -Xms20g"  # Java堆内存限制

# 线程和并行配置
export OMP_NUM_THREADS=4             # OpenMP线程数
export MKL_NUM_THREADS=4             # MKL线程数
export NUMEXPR_NUM_THREADS=4         # Numexpr线程数
export OPENBLAS_NUM_THREADS=4        # OpenBLAS线程数
export PYTHON_WORKER_THREADS=4       # Python工作线程
export SCENIC_NUM_WORKERS=10          # SCENIC工作进程数
export MAX_PARALLEL_SAMPLES=3        # 并行样本数

### —— 内存监控函数 —— ###
monitor_memory() {
    local mem_log="logs/memory_usage.log"
    echo "开始内存监控..." > "$mem_log"
    
    while true; do
        date >> "$mem_log"
        free -h >> "$mem_log"
        echo "进程内存使用:" >> "$mem_log"
        ps aux | grep -E 'python|R|java' | grep -v grep >> "$mem_log"
        echo "------------------------" >> "$mem_log"
        sleep 300  # 每5分钟记录一次
    done
}

### —— 资源检查函数 —— ###
check_resources() {
    local available_mem=$(free -g | awk '/^Mem:/{print $7}')
    local required_mem=160  # 需要至少160GB可用内存
    
    echo "系统资源检查:"
    echo "- 可用内存: ${available_mem}GB"
    echo "- 所需内存: ${required_mem}GB"
    echo "- 配置的最大使用: ${MAX_MEMORY_USAGE}GB"
    
    if [ $available_mem -lt $required_mem ]; then
        echo "警告: 可用内存不足"
        return 1
    fi
    return 0
}

### —— 目录结构初始化 —— ###
mkdir -p GRN_output/{raw_grn,final_results,stats_report}
mkdir -p logs/downstream

## —— 主程序 —— ###
main() {
    # 启动内存监控
    monitor_memory &
    MONITOR_PID=$!

    # 资源检查
    if ! check_resources; then
        echo "资源检查失败，退出运行"
        kill $MONITOR_PID
        exit 1
    fi

    # Stage 1: GRN运算
    echo "开始GRN运算..."
    for ((run=1; run<=NUM_RUNS; run++)); do
        echo "运行GRN迭代 ${run}"
        SEED=$((777 + run))

        arboreto_with_multiprocessing.py \
            input/mc_mat_for_step1.loom \
            ./cisTarget_db_new/mgi_hgnc_tfs.txt \
            --method grnboost2 \
            --output GRN_output/raw_grn/run_${run}_s1_adj.tsv \
            --num_workers $SCENIC_NUM_WORKERS \
            --seed $SEED > logs/grn_run_${run}.log 2>&1

        # 检查内存使用
        current_mem=$(free -g | awk '/^Mem:/{print $3}')
        if [ $current_mem -gt $MAX_MEMORY_USAGE ]; then
            echo "警告: 内存使用超过限制，暂停5分钟"
            sleep 300
        fi
    done

    # Stage 2: 运行合并脚本
    echo "运行网络合并脚本..."
    python 01_merge_script.py \
        --num-runs $NUM_RUNS \
        --importance-thresholds "$IMPORTANCE_THRESHOLDS" \
        --occurrence-thresholds "$OCCURRENCE_THRESHOLDS"

    # 清理监控进程
    kill $MONITOR_PID
}

# 运行主程序
main 2>&1 | tee logs/full_run_$(date +%Y%m%d_%H%M%S).log

### >>> Stage 2. Cistarget 运算 >>> ###
#!/bin/bash

### —— 设置基础参数 —— ###
# 设置基础路径 - 使用当前GRN输出目录
base_path="GRN_output/final_results"

# 在运行命令前添加环境变量
export PYTHONHASHSEED=42
export OMP_NUM_THREADS=4
export HDF5_USE_FILE_LOCKING=FALSE
export TMPDIR=`pwd`

# 激活conda环境（如果需要）
if [ -n "$CONDA_DEFAULT_ENV" ]; then
echo "Current conda environment: $CONDA_DEFAULT_ENV"
else
  echo "Activating conda environment..."
source ~/miniconda3/etc/profile.d/conda.sh
conda activate pyscenic_env
fi

# 检查Python依赖
check_python_dependencies() {
  echo "Checking Python dependencies..."
  python3 -c "
import sys
required_packages = ['pandas', 'loompy']
missing_packages = []
for package in required_packages:
    try:
        __import__(package)
        print(f'Found {package}')
    except ImportError:
        missing_packages.append(package)
if missing_packages:
    print('Missing packages: ' + ', '.join(missing_packages), file=sys.stderr)
    sys.exit(1)
"
  if [ $? -ne 0 ]; then
  echo "Error: Missing Python dependencies. Please install them first."
  exit 1
  fi
}

# 检查Python依赖
check_python_dependencies

# 固定输入文件
f_loom_grn="input/mc_mat_for_step1.loom"
f_loom_aucell="input/raw_count_mat_for_step3.loom"
f_tfs="./cisTarget_db_new/mgi_hgnc_tfs.txt"
f_motif_path="./cisTarget_db_new/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl"

# 创建日志目录
mkdir -p logs/downstream

# 检查必要文件是否存在的函数
check_files_exist() {
  local files=("$@")
  for file in "${files[@]}"; do
  if [ ! -f "$file" ]; then
  echo "Error: File $file does not exist"
  return 1
  fi
  done
  return 0
}

# 获取数据库文件列表的函数
get_database_files() {
  local db_path="./cisTarget_db_new"
  local db_files=(
    "$db_path/mm9-500bp-upstream-10species.mc9nr.genes_vs_motifs.rankings.feather"
    "$db_path/mm9-tss-centered-5kb-10species.mc9nr.genes_vs_motifs.rankings.feather"
    "$db_path/mm9-tss-centered-10kb-10species.mc9nr.genes_vs_motifs.rankings.feather"
  )
  
  # 检查文件是否存在
  check_files_exist "${db_files[@]}"
  
  # 如果文件存在，返回文件列表
  if [ $? -eq 0 ]; then
  printf '%s\n' "${db_files[@]}"
  else
    return 1
  fi
}

# 检查输入文件
echo "Checking input files..."
check_files_exist "$f_loom_grn" "$f_loom_aucell" "$f_tfs" "$f_motif_path"
if [ $? -ne 0 ]; then
echo "Error: Required input files are missing"
exit 1
fi

# 获取数据库文件列表
echo "Getting database files..."
f_db_names=$(get_database_files)
if [ $? -ne 0 ]; then
echo "Error: Database files not found"
exit 1
fi

# 设置随机种子
SEED=777

# 检查处理状态的函数
check_processing_status() {
  local param_dir="$1"
  local param_name=$(basename "$param_dir")
  local output_dir="$base_path/$param_name"
  
  # 检查所有必需的输出文件
  local required_files=(
    "${output_dir}/${param_name}_s2_reg.tsv"
    "${output_dir}/${param_name}_s3_aucell.loom"
    "${output_dir}/${param_name}_tf_stats.txt"
    "${output_dir}/${param_name}_s4_post_scenic/regulons.txt"
  )
  
  for file in "${required_files[@]}"; do
  if [ ! -f "$file" ]; then
  echo "Missing: $file"
  return 1
  fi
  done
  
  # 检查文件大小
  for file in "${required_files[@]}"; do
  if [ ! -s "$file" ]; then
  echo "Empty file: $file"
  return 1
  fi
  done
  
  return 0
}

# 修改process_network_file函数
process_network_file() {
  local network_file="$1"
  local param_dir=$(dirname "$network_file")
  local param_name=$(basename "$param_dir")
  
  # 创建输出目录
  local output_dir="$base_path/$param_name"
  mkdir -p "$output_dir"
  
  # 定义输出文件
  local ctx_output="$output_dir/${param_name}_s2_reg.tsv"
  local aucell_output="$output_dir/${param_name}_s3_aucell.loom"
  local post_scenic_output="$output_dir/${param_name}_s4_post_scenic"
  local stats_output="$output_dir/${param_name}_tf_stats.txt"
  
  # 创建日志文件
  local log_file="logs/downstream/${param_name}.log"
  
  # 检查处理状态
  if check_processing_status "$param_dir"; then
  echo "Skipping $param_name - all files complete and valid"
  return 0
  fi
  
  echo "Starting/Resuming processing for $param_name at $(date)" | tee -a "$log_file"
  
  # 检查并执行cisTarget
  if [ ! -f "$ctx_output" ] || [ ! -s "$ctx_output" ]; then
  echo "Running cisTarget..." | tee -a "$log_file"
  read -ra DB_FILES <<< "$f_db_names"
  
  pyscenic ctx \
  "$network_file" \
  "${DB_FILES[@]}" \
  --annotations_fname "$f_motif_path" \
  --expression_mtx_fname "$f_loom_grn" \
  --output "$ctx_output" \
  --num_workers 2 2>&1 | tee -a "$log_file"
  
  if [ $? -ne 0 ]; then
  echo "Error: cisTarget failed for $param_name" | tee -a "$log_file"
  return 1
  fi
  else
    echo "cisTarget output exists, skipping..." | tee -a "$log_file"
  fi
  
  # 检查并执行AUCell
  if [ ! -f "$aucell_output" ] || [ ! -s "$aucell_output" ]; then
  echo "Running AUCell..." | tee -a "$log_file"
  pyscenic aucell \
  "$f_loom_aucell" \
  "$ctx_output" \
  --output "$aucell_output" \
  --num_workers 2 \
  --seed "$SEED" 2>&1 | tee -a "$log_file"
  
  if [ $? -ne 0 ]; then
  echo "Error: AUCell failed for $param_name" | tee -a "$log_file"
  return 1
  fi
  else
    echo "AUCell output exists, skipping..." | tee -a "$log_file"
  fi
  
  # 更新统计信息
  echo "Updating statistics..." | tee -a "$log_file"
  python3 -c "
import loompy
import pandas as pd

# Initial TFs
initial_tfs = len(pd.read_csv('${network_file}', sep='\t')['TF'].unique())
print(f'Initial network TFs: {initial_tfs}')

# cisTarget TFs
cistarget_tfs = len(pd.read_csv('${ctx_output}', sep='\t')['TF'].unique())
print(f'After cisTarget TFs: {cistarget_tfs}')

# AUCell TFs
with loompy.connect('${aucell_output}') as ds:
    regulons = ds.ra['RegulonName']
    unique_tfs = set([reg.split('_')[0] for reg in regulons])
    print(f'After AUCell TFs: {len(unique_tfs)}')
" > "$stats_output" 2>> "$log_file"
  
  # 检查并执行post-processing
  if [ ! -f "${post_scenic_output}/regulons.txt" ] || [ ! -s "${post_scenic_output}/regulons.txt" ]; then
  echo "Running post-processing..." | tee -a "$log_file"
  local threads=1
  local min_regulon_size=10
  
  bash 02_postSCENIC.sh "$aucell_output" "$ctx_output" "$post_scenic_output" "$threads" "$min_regulon_size" 2>&1 | tee -a "$log_file"
  
  if [ $? -ne 0 ]; then
  echo "Error: Post-SCENIC processing failed for $param_name" | tee -a "$log_file"
  return 1
  fi
  
  # 添加最终regulon统计
  if [ -f "${post_scenic_output}/regulons.txt" ]; then
  final_tf_count=$(cut -f1 "${post_scenic_output}/regulons.txt" | sort -u | wc -l)
  echo "Final regulon TFs: $final_tf_count" >> "$stats_output"
  fi
  else
    echo "Post-processing output exists, skipping..." | tee -a "$log_file"
  fi
  
  echo "Completed processing for $param_name at $(date)" | tee -a "$log_file"
  return 0
}

# 导出函数和变量，使其在子进程中可用
export -f process_network_file check_files_exist
export f_db_names f_loom_grn f_loom_aucell f_motif_path SEED
export base_path

# 定义并行数量
max_parallel_samples=4

# 存储后台进程PID
pids=()

# 查找所有 network.tsv 文件
echo "Finding network files..."
network_files=($(find "$base_path" -name "network.tsv"))
total_files=${#network_files[@]}
  echo "Found $total_files network files to process"
  
  # 创建主日志文件
  exec 1> >(tee -a "logs/downstream/main.log") 2>&1
  
  # 并行处理网络文件
  echo "Starting parallel processing of network files..."
  for ((i=0; i<total_files; i+=max_parallel_samples)); do
  echo "Processing batch starting at index $i"
  pids=()
  
  # 处理当前批次的文件
  for ((j=0; j<max_parallel_samples && i+j<total_files; j++)); do
  process_network_file "${network_files[i+j]}" &
    pid=$!
    pids+=($pid)
  echo "Started process $pid for ${network_files[i+j]}"
  done
  
  # 等待当前批次完成
  for pid in "${pids[@]}"; do
  wait "$pid"
  status=$?
    if [ $status -ne 0 ]; then
  echo "Process $pid failed with status $status"
  else
    echo "Process $pid completed successfully"
  fi
  done
  
  echo "Completed batch processing from index $i"
  done
  
  echo "All network processing completed at $(date)"
  
  # 生成处理报告
  echo "Generating processing report..."
  echo "Processing Summary" > "logs/downstream/processing_report.txt"
  echo "==================" >> "logs/downstream/processing_report.txt"
  echo "Total files processed: $total_files" >> "logs/downstream/processing_report.txt"
  echo "Timestamp: $(date)" >> "logs/downstream/processing_report.txt"
  echo "" >> "logs/downstream/processing_report.txt"
  
  # 检查每个输出目录的状态并添加TF统计信息
  for network_file in "${network_files[@]}"; do
  param_dir=$(dirname "$network_file")
  param_name=$(basename "$param_dir")
  output_dir="$base_path/$param_name"
  stats_file="$output_dir/${param_name}_tf_stats.txt"
  
  echo "Parameter set: $param_name" >> "logs/downstream/processing_report.txt"
  if [ -f "$output_dir/${param_name}_s3_aucell.loom" ]; then
  echo "Status: Complete" >> "logs/downstream/processing_report.txt"
  if [ -f "$stats_file" ]; then
  echo "TF Statistics:" >> "logs/downstream/processing_report.txt"
  cat "$stats_file" >> "logs/downstream/processing_report.txt"
  fi
  else
    echo "Status: Incomplete or Failed" >> "logs/downstream/processing_report.txt"
  fi
  echo "-------------------" >> "logs/downstream/processing_report.txt"
  done
  
  # 生成TF统计总结
  echo "" >> "logs/downstream/processing_report.txt"
  echo "TF Statistics Summary" >> "logs/downstream/processing_report.txt"
  echo "====================" >> "logs/downstream/processing_report.txt"
  echo "Parameter Set | Initial TFs | After cisTarget | After AUCell | Final Regulon" >> "logs/downstream/processing_report.txt"
  echo "-------------|-------------|-----------------|--------------|---------------" >> "logs/downstream/processing_report.txt"
  
  for network_file in "${network_files[@]}"; do
  param_dir=$(dirname "$network_file")
  param_name=$(basename "$param_dir")
  stats_file="$base_path/$param_name/${param_name}_tf_stats.txt"
  
  if [ -f "$stats_file" ]; then
  initial_tfs=$(grep "Initial network TFs:" "$stats_file" | cut -d' ' -f4)
  cistarget_tfs=$(grep "After cisTarget TFs:" "$stats_file" | cut -d' ' -f4)
  aucell_tfs=$(grep "After AUCell TFs:" "$stats_file" | cut -d' ' -f4)
  final_tfs=$(grep "Final regulon TFs:" "$stats_file" | cut -d' ' -f4)
  
  echo "$param_name | $initial_tfs | $cistarget_tfs | $aucell_tfs | $final_tfs" >> "logs/downstream/processing_report.txt"
  fi
  done
  
  echo "" >> "logs/downstream/processing_report.txt"
  echo "Processing report generated at logs/downstream/processing_report.txt" 