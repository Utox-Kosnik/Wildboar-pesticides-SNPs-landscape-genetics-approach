library(dplyr)
library(tidyr)
library(purrr)
library(data.table)

dg_getHypergeomAssociation_v6 <- function(df1, df2, colname1, colname2, p_adj_method = "fdr", background = NULL, showParam = FALSE) {
  
  # ids or unique identifiers need to be in the first column of df1 and df2
  # check if df1 and df2 are in the right format
  # check if sensitivity value is a vector
  
  # df1 <- ctd_gene_chem_interest
  # df2 <- ctd_gene_disease
  # 
  # colname1 <- "GeneSymbol"
  # colname2 <- "X..GeneSymbol"
  # overlap_sensitivity <- c(3)
  # p_adj_method <- "bonferroni"
  
  # df1 <- gene_phenotype
  # df2 <- gene_pathway
  # colname1 <- "gene"
  # colname2 <- "gene"
  # p_adj_method <- "bonferroni"
  # background <- NULL
  # 
  # df1 <- path.dis 
  # df2 <- chem.gene
  # colname1 <- "GeneSymbol"
  # colname2 <- "GeneSymbol"
  # overlap_sensitivity <- c(1)
  # p_adj_method <- "fdr"
  
  # df1 <- gene_phenotype
  # df2 <- gene_pathway
  # colname1 <- "gene"
  # colname2 <- "gene"
  # p_adj_method <- "fdr"
  
  setDT(df1)
  setDT(df2)
  
  coldf1 <- names(df1)[1]
  coldf2 <- names(df2)[1]
  
  
  df1_list <- df1 %>%
    group_by(!!sym(coldf1)) %>%
    summarize(!!paste0("genes_", coldf1) := list(unique(!!sym(colname1)))) %>%
    ungroup()
  df2_list <- df2 %>%
    group_by(!!sym(coldf2)) %>%
    summarize(!!paste0("genes_", coldf2) := list(unique(!!sym(colname2)))) %>%
    ungroup()
  
  # creating unique combination of df1_ids and df2_ids
  cat("Creation of unique combination datatable \n")
  combinations_dt <- CJ(unique(df1[[1]]), unique(df2[[1]]))
  setnames(combinations_dt, c(colnames(df1)[1], colnames(df2)[1]))
  combinations_dt[, unique_id := paste(get(coldf1), get(coldf2), sep = "_")]
  
  #qvalues overlap genes
  cat("qvalues calculation \n")
  q.dataset <- merge(df1[, .(gene = get(colname1), c1 = get(coldf1))], 
                     df2[, .(gene = get(colname2), c2 = get(coldf2))], 
                     by = "gene", allow.cartesian = TRUE)
  q.dataset$unique_id <- paste0(q.dataset$c1, "_", q.dataset$c2)
  q.dataset<- q.dataset[, .(num_gene_overlap = length(unique(gene))), by=.(unique_id)] #count the unique gene per unique_id (unique pathway_phenotype combination) hence this is gene overlap
  q_all<- q.dataset[match(combinations_dt$unique_id, q.dataset$unique_id)]
  
  combinations_dt$qvalue <- q_all$num_gene_overlap
  combinations_dt[is.na(qvalue)]$qvalue<- 0
  
  combinations_dt$qvalue <- combinations_dt$qvalue - 1 ### FIX, we substract -1 to the qvalue in the hypergeom test
  
  # here could filter out 0 gene overlap
  combinations_dt <- combinations_dt[combinations_dt$qvalue > 0,]
  
  #mvalues: number of unique genes associated with each df1 code
  cat("mvalues calculation \n")
  m.vals <- df1[, .(genes = length(unique(get(colname1)))), by = .(get(coldf1))]
  m_all <- m.vals[match(combinations_dt[[1]], m.vals[[1]])]
  combinations_dt$mvalue <- m_all$genes
  
  #nvalues:  number of genes total – number of genes in a df1
  cat("nvalue calculation \n")
  if (is.null(background)) {
    cat("BEWARE! The background used is the combination of elements between the 2 dataset, this can depends on your question \n")
    
    background <- length(unique(c(df1[[colname1]], df2[[colname2]])))
  }
  combinations_dt$nvalue <- background - combinations_dt$mvalue
  
  #kvalues: number of unique genes associated with each pathway code
  cat("kvalue calculation \n")
  k.vals <- df2[, .(genes = length(unique(get(colname2)))), by = .(get(coldf2))]
  k_all <- k.vals[match(combinations_dt[[2]], k.vals[[1]])]
  combinations_dt$kvalue <- k_all$genes
  
  #pvalue calculation
  cat("pvalue calculation \n")
  combinations_dt$pvalue <- phyper(q=combinations_dt$qvalue, m=combinations_dt$mvalue, n=combinations_dt$nvalue, k=combinations_dt$kvalue, lower.tail=FALSE)
  combinations_dt$p_value <- p.adjust(combinations_dt$pvalue, method = p_adj_method) # adjust the pvalue
  combinations_dt <- combinations_dt[combinations_dt$p_value < 0.05, ]
  
  cat("Adding final informations to the association datatable \n")
  association <- as.data.table(combinations_dt)
  
  if (!showParam) {
    association[, c("qvalue", "mvalue", "nvalue", "kvalue", "pvalue") := NULL]
  }
  
  association <- merge(association, df1_list, by = coldf1, all.x = TRUE)
  association <- merge(association, df2_list, by = coldf2, all.x = TRUE)
  association <- merge(association, q_all[, .(unique_id, num_gene_overlap)], by = "unique_id", all.x = TRUE)
  
  
  return(association)
}
