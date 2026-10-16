############################################################################################
################################# RA + GENE ASSIGNMENT #####################################
############################################################################################

Sys.setenv(LANG = "en")
rm(list = ls()) 

######################
#Packages needed (some packages are called within the code)
######################

library(dplyr) #%>%
#other packages specified in code: vegan, qvalue, rtracklayer
#packageVersion("rtracklayer")
#citation("rtracklayer")
#packageVersion("vegan")
#citation("vegan")
#packageVersion("qvalue")
#citation("qvalue")


######################
#Input files needed
######################

#CHEMICAL_LIST
chemical_list <- readRDS("chemical_list.RData") #FROM SCRIPT 2

#MATRIX FILTERED FOR 282 INDIVIDUALS
wb_matrix_no_na_282 <- readRDS("wb_matrix_no_na_282.RData") #FROM SCRIPT 2

#FINAL METADATA
wb_final_metadata <-readRDS("wb_final_metadata_4km.RData") #FROM SCRIPT 2
str(wb_final_metadata)

#ILLUMINA DOCUMENT:match SNP name with real chromosome + position (to later version of the reference genome: from 10.2 to 11.1)
snp_id<- read.table("forR_PorcineSNP60.csv",  sep=",", header = T) #FROM ILLUMINA

#SNP names, taken from map file
all_snp_names_from_map <- read.table("names_SNPs.txt", header = F) #FROM MAP file

#GENOME DATA TO MATCH SNPS WITH GENES
g <- rtracklayer::readGFF("Ensembl_Sus_scrofa.Sscrofa11.1.112.gtf") #FROM ENSEMBL
str(g)
g$start<-as.numeric(g$start)
g$end <- as.numeric(g$end)

test_match_order2 <- function(x,y) {
  if (isTRUE(all.equal(x,y))) print('Perfect match in same order')
  if (!isTRUE(all.equal(x,y)) && isTRUE(all.equal(sort(x),sort(y)))) print('Perfect match in wrong order')
  if (!isTRUE(all.equal(x,y)) && !isTRUE(all.equal(sort(x),sort(y)))) print('No match')
}

test_match_order2(rownames(wb_matrix_no_na_282), (wb_final_metadata$id)) #"Perfect match in same order"

############################################################################################

######################
#Run the complete for-loop
######################

all_tox_wb = NULL
tox_without_snps = NULL
all_vif_results <- NULL
all_rda_results <- NULL

for (chemical in chemical_list){
  print(chemical)
  toxic_column <- unlist(wb_final_metadata[chemical]) 
  
  set.seed(123)
  wb.rda.tox <- vegan::rda(wb_matrix_no_na_282 ~ toxic_column + Condition(longitude+latitude+agricultural_percentage+human_population), data=wb_final_metadata, scale = T, na.action = na.omit) 
  
  ############################################################################################
  
  ###check overinflation of the model (VIF)
  
  vif_result <- tryCatch(
    vegan::vif.cca(wb.rda.tox),
    error = function(e) {
      message("Skipping ", chemical, ": pesticide is redundant with the conditional variables.", conditionMessage(e))
      return(NULL)
    }
  )
  
  if (is.null(vif_result)) next
  
  print(vif_result)
  
  #ONLY 2 CHEMICALS ARE PROBLEMATIC
  #Bifenthrin 
  #Forchlorfenuron
  
  # store VIF results
  vif_df <- data.frame(
    chemical = chemical,
    variable = names(vif_result),
    VIF = as.numeric(vif_result)
  )
  
  all_vif_results <- rbind(all_vif_results, vif_df)
  
  ############################################################################################
  ### Overall significance and explained variation of pesticide
  
  # permutation test
  set.seed(123)
  rda_test <- anova(wb.rda.tox, permutations = 999)
  
  # R2
  r2_result <- vegan::RsquareAdj(wb.rda.tox)
  
  # store results
  rda_df <- data.frame(
    chemical = chemical,
    F_value = rda_test$F[1],
    p_value = rda_test$`Pr(>F)`[1],
    R2 = r2_result$r.squared,
    adjusted_R2 = r2_result$adj.r.squared,
    variation_explained_percent = r2_result$r.squared * 100,
    adjusted_variation_explained_percent = r2_result$adj.r.squared * 100
  )
  
  all_rda_results <- rbind(all_rda_results, rda_df)
  
  ############################################################################################
  
  ### p-values
  
  load.rda.class <- vegan::scores(wb.rda.tox, choices=c(1), display="species") #species = SNPs, sites = individuals
  load.rda.class = as.data.frame(load.rda.class)
  load.rda.class$locus = rownames(load.rda.class)
  #group loadings from all axes together
  mload.rda.class = tidyr::gather(load.rda.class, axis, loading, -locus )
  #mean-center and scale loadings to variance = 1
  mload.rda.class$scaled_loading = scale(mload.rda.class$loading) #zscore range goes from -4 until 4
  mload.rda.class$pvalues <- (1-pnorm(abs(mload.rda.class$scaled_loading)))*2 #zscore to two-tailed p-value 
  
  ### Calculate q-values
  
  qvalues.class <- qvalue::qvalue(mload.rda.class$pvalues)
  mload.rda.class$qvalues <- qvalues.class$qvalues
  df_rda.qvalues.class <- as.data.frame(mload.rda.class)
  snp_wb_tox <- subset(df_rda.qvalues.class, qvalues < 0.05) 
  
  ######################################################################################################
  
  #Significant SNPs
  loci_wb_tox <- snp_wb_tox$locus #
  
  df_rda_final <- snp_wb_tox[,c(1,6)]
  colnames(df_rda_final)<-c("snp","qvalue") 
  
  ########################################################################################################
  #MATCH SNPS TO CHROMOSOME + POSITION: WILD BOAR
  #######################################################################################################
  

  list_of_values <- df_rda_final$snp
  list_of_qvalues <- df_rda_final$qvalue
  all_snp_names_from_map$number <- rownames(all_snp_names_from_map )
  colnames(all_snp_names_from_map )<-c("chrom","name","x","pos","nr") #this comes from the map file --> BUT chrom + pos are wrong since this comes from a different reference genome!
  filtered_df <- subset(all_snp_names_from_map , nr %in% list_of_values)
  snp_names <- filtered_df$name
  
  all_snp = NULL
  i=1
  for (snp in snp_names){
    help <- snp_id %>% dplyr::filter(snp_id$Name == snp)%>% as.data.frame()
    chromosome <- help$Chr
    position <- help$MapInfo
    q <- list_of_qvalues[i]
    new_snp <- c(snp,chromosome,position, q)
    all_snp <- rbind (all_snp, new_snp)
    i = i + 1
  }
  
  all_snp <-as.data.frame(all_snp)
  
  ########################################################################################################
  #MATCH SNPS TO GENES
  ####################################################################################################### 
  
  if (nrow(all_snp) > 0 ){
    colnames(all_snp) <- c("snp","chromosome","position","qvalue")
    all_snp$position <- as.numeric(all_snp$position)
    rownames(all_snp) <- c(1:(nrow(all_snp)))  
    n = 0
    i=1
    final_df = NULL
    for(i in 1:nrow(all_snp)){
      print(i)
      right_snp <- all_snp$snp[i]
      right_chromosome <- all_snp$chr[i]
      right_position <- all_snp$position[i]
      qvalue <- all_snp$qvalue[i]
      g_filter <- g %>% dplyr::filter(seqid == right_chromosome, type == "gene", 
                               (start <= right_position & right_position <= end)) %>% as.data.frame()  #intragenic SNPs
      g_filter <- g_filter[!duplicated(g_filter), ]
      if (nrow(g_filter) == 0){ #no closeby genes
        g_subset <- g_filter[nrow(g_filter)+1,] <- NA
        g_subset <- g_filter[,c(1:5,9,11)]
      } else { #intragenic genes
        n <- n + 1
        g_subset <-g_filter[,c(1:5,9,11)]
      }
      g_subset$snp <- right_snp
      g_subset$chromosome <- right_chromosome
      g_subset$position <- right_position
      g_subset <- cbind(g_subset, qvalue)
      final_df = rbind(final_df,g_subset)
      
      i = i + 1
    }
    reduced_df <- final_df[,c(8:10,7,6,4,5,11)] #0 kB away
    reduced_df$tox_name <- chemical
    all_tox_wb <- rbind(all_tox_wb,reduced_df) 
  } else {
    tox_without_snps <- cbind(tox_without_snps, chemical) #
  }}



unique(all_tox_wb$tox_name) #

saveRDS(all_tox_wb, "all_tox_wb_ref_11_1_revised_4km.RData") 

###############################

unique(all_vif_results$chemical) #120 (Bifenthrin was excluded)
all_vif_results$VIF <- as.integer(round(all_vif_results$VIF, 0)) #should be fine (9.6 can be turned into 10 but checked results and we are good)
all_vif_results_wide <- all_vif_results %>%tidyr::pivot_wider(names_from = variable,values_from = VIF)
final_model_results <- all_vif_results_wide %>%left_join(all_rda_results, by = "chemical")

write.csv(final_model_results,"vif_results_climate_excluded_4km.csv",row.names = FALSE)

final_model_results_summary <- final_model_results %>%
  mutate(class = case_when(
      if_any(c(longitude, latitude, agricultural_percentage,human_population, toxic_column),~ .x >= 10) ~ "VIF >= 10",
      if_all(c(longitude, latitude, agricultural_percentage,human_population, toxic_column),~ .x < 10) & p_value > 0.05 ~ "VIF < 10 + p > 0.05",
      if_all(c(longitude, latitude, agricultural_percentage,human_population, toxic_column),~ .x < 10) & p_value <= 0.05 ~ "VIF < 10 + p <= 0.05"))
table(final_model_results_summary$class)

# Chemicals with acceptable VIF and significant overall RDA
significant_model_chemicals <- final_model_results_summary %>%
  filter(class == "VIF < 10 + p <= 0.05") %>%
  pull(chemical) %>%
  unique()

# Chemicals with significant SNPs
significant_snp_chemicals <- unique(all_tox_wb$tox_name)

# Are the two sets exactly identical?
setequal(significant_model_chemicals, significant_snp_chemicals)

# Present in all_tox_wb, but overall model is not in the significant class
setdiff(significant_snp_chemicals,significant_model_chemicals)
#9 chemicals have SNPs but no significant overall model

