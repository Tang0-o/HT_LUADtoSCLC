#!/bin/bash

# 接收参数
aucell_output=$1
ctx_output=$2
sample_name=$3
threads=$4
min_regulon_size=$5

# 执行脚本
python 02_postSCENIC.py \
$aucell_output \
$ctx_output \
$sample_name \
$threads \
$min_regulon_size