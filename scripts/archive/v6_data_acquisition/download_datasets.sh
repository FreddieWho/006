#!/bin/bash
# 006项目数据下载脚本
# 生成时间: 2026-05-05

set -e

BASE_DIR="/home/huyudi/006/data/raw"
LOG_FILE="/home/huyudi/006/docs/download_log.txt"

echo "=== 006项目数据下载开始 ===" | tee -a "$LOG_FILE"
echo "开始时间: $(date)" | tee -a "$LOG_FILE"

# GEO数据下载函数
download_geo() {
    local gse_id=$1
    local data_type=$2
    local target_dir="$BASE_DIR/$data_type/GEO/$gse_id"
    
    echo "下载 $gse_id ($data_type)..." | tee -a "$LOG_FILE"
    mkdir -p "$target_dir"
    
    # 获取GSE编号的前缀用于FTP路径
    local gse_prefix=$(echo $gse_id | sed 's/GSE\([0-9]*\).*/\1/' | sed 's/\(.*\)\(nnn\)/\1/')
    local gse_series="GSE${gse_prefix:0:3}nnn"
    
    # 下载
    cd "$target_dir"
    wget -r -np -nH --cut-dirs=5 \
        "ftp://ftp.ncbi.nlm.nih.gov/geo/series/$gse_series/$gse_id/" \
        2>&1 | tee -a "$LOG_FILE" || echo "警告: $gse_id 下载可能不完整" | tee -a "$LOG_FILE"
    
    echo "$gse_id 下载完成" | tee -a "$LOG_FILE"
}

# SRA数据下载函数
download_sra() {
    local bioproject_id=$1
    local target_dir="$BASE_DIR/scRNA/SRA/$bioproject_id"
    
    echo "下载 $bioproject_id (SRA)..." | tee -a "$LOG_FILE"
    mkdir -p "$target_dir"
    cd "$target_dir"
    
    # 使用prefetch下载
    prefetch "$bioproject_id" 2>&1 | tee -a "$LOG_FILE" || \
        echo "警告: $bioproject_id prefetch失败，尝试直接下载" | tee -a "$LOG_FILE"
    
    echo "$bioproject_id 下载完成" | tee -a "$LOG_FILE"
}

# 1. GSE236581 - 结直肠癌新辅助抗PD-1 (scRNA+scTCR)
download_geo "GSE236581" "scRNA"

# 2. GSE221561 - 食管鳞癌PD-1+化疗/放疗 (scRNA)
download_geo "GSE221561" "scRNA"

# 3. GSE256326 - 肾髓质癌nivolumab+ipilimumab (scRNA)
download_geo "GSE256326" "scRNA"

# 4. GSE314072 - 肾细胞癌 (scRNA+scTCR, 参考)
download_geo "GSE314072" "scRNA"

# 5. GSE289745 - 皮肤鳞癌Visium (ST)
download_geo "GSE289745" "ST"

# 6. GSE291246 - 基底细胞癌Xenium (ST)
download_geo "GSE291246" "ST"

# 7. GSE235672 - 胶质母细胞瘤Visium (ST)
download_geo "GSE235672" "ST"

# 8. GSE238264 - 肝癌Visium (ST)
download_geo "GSE238264" "ST"

# 9. GSE177043 - 三阴性乳腺癌Visium (ST)
download_geo "GSE177043" "ST"

# 10. PRJNA932556 - 结直肠癌MSI-H抗PD-1 (scRNA, SRA)
download_sra "PRJNA932556"

echo "" | tee -a "$LOG_FILE"
echo "=== 下载完成 ===" | tee -a "$LOG_FILE"
echo "完成时间: $(date)" | tee -a "$LOG_FILE"
echo "注意: OMIX005710和PRJCA039752需要手动从网站下载" | tee -a "$LOG_FILE"
echo "注意: HRA003591, HRA007492, HRA007492_ST需要申请Controlled Access" | tee -a "$LOG_FILE"
echo "注意: NCT02451982需要从原文获取具体GEO/SRA编号" | tee -a "$LOG_FILE"
