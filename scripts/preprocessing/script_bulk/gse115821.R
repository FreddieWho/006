
#  loading cts ------------------------------------------------------------

cts <- data.table::fread('~/006/data/combo/GSE115821/GSE115821_MGH_counts.csv',
                         check.names = F)


library(xml2)
meta_xml <- read_xml('~/006/data/combo/GSE115821/GSE115821_family.xml')
meta <- data.frame(
  id = xml_attr(xml_find_all(xml_data, "//book"), "id"),
  title = xml_text(xml_find_all(xml_data, "//title")),
  author = xml_text(xml_find_all(xml_data, "//author")),
  price = as.numeric(xml_text(xml_find_all(xml_data, "//price")))
)

ns <- xml_ns(meta_xml)

extract_samples <- function(xml_data) {
  samples <- xml_find_all(xml_data, "//d1:Sample", ns)
  
  result <- lapply(samples, function(s) {
    # 基本信息
    iid <- xml_attr(s, "iid")
    accession <- xml_text(xml_find_first(s, ".//d1:Accession", ns))
    title <- xml_text(xml_find_first(s, ".//d1:Title", ns))
    type <- xml_text(xml_find_first(s, ".//d1:Type", ns))
    
    # Channel 信息
    source <- xml_text(xml_find_first(s, ".//d1:Source", ns))
    organism <- xml_text(xml_find_first(s, ".//d1:Organism", ns))
    molecule <- xml_text(xml_find_first(s, ".//d1:Molecule", ns))
    
    # Platform
    platform_ref <- xml_attr(xml_find_first(s, ".//d1:Platform-Ref", ns), "ref")
    
    # Library 信息
    library_strategy <- xml_text(xml_find_first(s, ".//d1:Library-Strategy", ns))
    library_source <- xml_text(xml_find_first(s, ".//d1:Library-Source", ns))
    instrument <- xml_text(xml_find_first(s, ".//d1:Predefined", ns))
    
    # Characteristics (特征信息)
    chars <- xml_find_all(s, ".//d1:Characteristics", ns)
    char_tags <- xml_attr(chars, "tag")
    char_values <- trimws(xml_text(chars))
    
    # 创建基础数据框
    base_df <- data.frame(
      Sample_ID = iid,
      Accession = accession,
      Title = title,
      Type = type,
      Source = source,
      Organism = organism,
      Molecule = molecule,
      Platform = platform_ref,
      Library_Strategy = library_strategy,
      Library_Source = library_source,
      Instrument = instrument,
      stringsAsFactors = FALSE
    )
    
    # 添加 Characteristics
    if (length(char_tags) > 0) {
      char_df <- as.data.frame(t(char_values))
      colnames(char_df) <- char_tags
      base_df <- cbind(base_df, char_df)
    }
    
    return(base_df)
  })
  
  do.call(bind_rows, result)
}

samples_df <- extract_samples(meta_xml)
rownames(samples_df) <- samples_df$Title

colnames(cts) <- gsub('\\.bam','',colnames(cts))

gse115821 <- list('exp' = cts,
                  'meta' = samples_df)

qs::qsave(gse115821,file = '~/006/data/processed/rds/bulkrna_gse115821.qs')
