library(rtracklayer)
library(dplyr)
library(tidyr)
library(Seurat)
# -------------------------
# Step 1: 读取离线Gencode gtf，提取gene_id, gene_name, gene_type
# read_gencode_gtf <- function(gtf_path) {
#   gtf <- import(gtf_path)
#   genes <- gtf[gtf$type == "gene"]
#   df <- data.frame(
#     gene_id = mcols(genes)$gene_id,
#     gene_name = mcols(genes)$gene_name,
#     gene_type = mcols(genes)$gene_type,
#     stringsAsFactors = FALSE
#   )
#   df <- df %>% filter(!is.na(gene_name) & gene_name != "") %>% distinct()
#   return(df)
# }
# 
# gencode_v19 <- read_gencode_gtf("~/006/data/ref/gencode.v49lift37.basic.annotation.gtf")
# gencode_v38 <- read_gencode_gtf("~/006/data/ref/gencode.v49.basic.annotation.gtf")
# gene_anno_list <- list(v19 = gencode_v19, v38 = gencode_v38)
# saveRDS(gene_anno_list,file = '~/006/data/ref/gene_anno_list.rds')

gene_anno_list <- readRDS('~/006/data/ref/gene_anno_list.rds')

# Step 2: 读取HGNC别名文件，构建alias->approved symbol 映射表
read_hgnc_alias <- function(hgnc_file) {
  hgnc <- read.delim('~/006/data/ref/hgnc_complete_set.txt', stringsAsFactors = FALSE, quote = "", check.names = FALSE)
  
  library(dplyr)
  library(tidyr)
  
  # 处理 alias_symbol 列，拆分成多行，忽略空
  alias1 <- hgnc %>%
    select(symbol, alias_symbol) %>%
    mutate(alias_symbol = ifelse(is.na(alias_symbol), "", alias_symbol)) %>%
    filter(alias_symbol != "") %>%
    mutate(alias_symbol = gsub('"','',alias_symbol)) %>% 
    separate_rows(alias_symbol, sep = "\\|") %>%
    rename(alias = alias_symbol) %>%
    mutate(alias = trimws(alias))
  
  # 处理 prev_symbol 列，拆分成多行，忽略空
  alias2 <- hgnc %>%
    select(symbol, prev_symbol) %>%
    mutate(prev_symbol = ifelse(is.na(prev_symbol), "", prev_symbol)) %>%
    filter(prev_symbol != "") %>%
    mutate(prev_symbol = gsub('"','',prev_symbol)) %>% 
    separate_rows(prev_symbol, sep = "\\|") %>%
    rename(alias = prev_symbol) %>%
    mutate(alias = trimws(alias))
  
  # 合并两个数据框，去重
  alias_df <- bind_rows(alias1, alias2) %>%
    distinct(symbol,alias)
  
  # 构建别名到官方symbol映射的named vector
  alias_map <- setNames(alias_df$symbol, alias_df$alias)
  
  return(alias_map)
}

alias_map <- read_hgnc_alias('~/006/data/ref/hgnc_complete_set.txt')
saveRDS(alias_map,file = '~/006/data/ref/alias_map.rds')
# Step 3: 识别基因ID类型
# 返回值为字符标识： "ensembl" / "entrez" / "symbol"
# 简单规则基于格式和数字比例
detect_id_type <- function(ids) {
  sapply(ids, function(id) {
    if (grepl("^ENS[A-Z]*G[0-9]+$", id)) {
      return("ensembl")
    } else if (grepl("^[0-9]+$", id)) {
      return("entrez")
    } else {
      # 基于是否为符号（一般含字母且非数字全部）
      return("symbol")
    }
  }, USE.NAMES = FALSE)
}

# Step 4: 构建不同ID转换映射表函数
# 输入：ids，注释表，alias表
# 优先使用注释表映射ids到symbol，再用alias映射转换别名
convert_to_symbol <- function(ids, id_type, gene_anno_list, alias_map) {
  gene_v19 <- gene_anno_list$v19
  gene_v38 <- gene_anno_list$v38
  
  symbol_out <- character(length(ids))
  
  # 分别处理不同类型的ids
  ens_idx <- which(id_type == "ensembl")
  ent_idx <- which(id_type == "entrez")
  sym_idx <- which(id_type == "symbol")
  
  # Ensembl ID优先v38，再v19
  if(length(ens_idx) > 0) {
    ens_ids <- ids[ens_idx]
    map38 <- gene_v38$gene_name[match(ens_ids, gene_v38$gene_id)]
    miss38 <- is.na(map38)
    if(any(miss38)) {
      map19 <- gene_v19$gene_name[match(ens_ids[miss38], gene_v19$gene_id)]
      map38[miss38] <- map19
    }
    symbol_out[ens_idx] <- map38
  }
  
  # Entrez ID映射留空或扩展（这里先NA）
  if(length(ent_idx) > 0) {
    symbol_out[ent_idx] <- NA_character_
  }
  
  # Symbol类型：先判断是否在官方symbol里（v38优先，下同）
  if(length(sym_idx) > 0) {
    sym_ids <- ids[sym_idx]
    in_v38 <- sym_ids %in% gene_v38$gene_name
    in_v19 <- sym_ids %in% gene_v19$gene_name
    alias_in <- sym_ids %in% names(alias_map)
    
    # 预留输出向量
    sym_res <- character(length(sym_ids))
    
    sym_res[in_v38] <- sym_ids[in_v38]
    sym_res[!in_v38 & in_v19] <- sym_ids[!in_v38 & in_v19]
    sym_res[!in_v38 & !in_v19 & alias_in] <- alias_map[sym_ids[!in_v38 & !in_v19 & alias_in]]
    
    # 剩下未匹配的保留原符号
    unmatched <- !(in_v38 | in_v19 | alias_in)
    sym_res[unmatched] <- sym_ids[unmatched]
    
    symbol_out[sym_idx] <- sym_res
  }
  
  return(symbol_out)
}
# Step 5: 根据symbol匹配gene_type（优先v38，其次v19）
get_gene_type <- function(symbols, gene_anno_list) {
  gene_v19 <- gene_anno_list$v19
  gene_v38 <- gene_anno_list$v38
  
  # 先在v38匹配
  idx38 <- match(symbols, gene_v38$gene_name)
  type38 <- gene_v38$gene_type[idx38]
  
  # 对v38未匹配的，转去v19匹配
  missing <- is.na(type38)
  idx19 <- match(symbols[missing], gene_v19$gene_name)
  type19 <- gene_v19$gene_type[idx19]
  
  # 合并结果
  type38[missing] <- type19
  
  return(type38)
}
# ------------------------

gene_symbol_align <- function(input_ids,
                              gene_anno_list,
                              alias_map) {

  message("Detecting input ID types...")
  id_types <- detect_id_type(input_ids)
  
  message("Mapping IDs to approved gene symbols...")
  aligned_symbols <- convert_to_symbol(input_ids, id_types, gene_anno_list, alias_map)
  
  message("Getting gene biotype annotations...")
  gene_types <- get_gene_type(aligned_symbols, gene_anno_list)
  
  # 结果表
  out_df <- data.frame(
    original_id = input_ids,
    aligned_symbol = aligned_symbols,
    gene_type = gene_types,
    conversion_status = ifelse(is.na(aligned_symbols), "failed", "success"),
    stringsAsFactors = FALSE
  )
  
  return(out_df)
}

# === 使用示例 ===
# 请替换路径为你的本地文件路径
gene_anno_list <- readRDS('~/006/data/ref/gene_anno_list.rds')
alias_map <- readRDS(file = '~/006/data/ref/alias_map.rds')
# 输入基因ID样例（可替换成你的输入基因向量）
# input_ids <- c("ENSG00000171554", "IGJ", "JCHAIN", "7157", "TP53", "XYZ")
