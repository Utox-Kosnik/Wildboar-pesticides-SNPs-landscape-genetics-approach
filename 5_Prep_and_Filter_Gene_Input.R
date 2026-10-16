############################################################################################
###################### PREP AND FILTER GENE INPUT #########################################
############################################################################################

Sys.setenv(LANG = "en")
rm(list = ls()) 

######################
#Packages needed (some packages are called within the code)
######################

library(dplyr) #%>%
citation("data.table")
packageVersion("data.table")
######################
#Input files needed
######################

#Output from RDA analysis
all_tox_wb <- readRDS("all_tox_wb_ref_11_1_revised_4km.RData") #FROM SCRIPT 3
unique(all_tox_wb$tox_name)

#looked up the protein linked to the ENSSSCG codes and merged it into a text file
ENS_names <- read.table("ENS_codes_with_names_0kb_11_1.txt", header = T) #FROM SCRIPT 3, modified out of R

#CHEMICAL_LIST
chemical_list <- readRDS("chemical_list.RData") #FROM SCRIPT 2

####################################################################################################################################################################
################################ADD ENS GENE NAMES####################################################################################################################################
####################################################################################################################################################################

all_tox_wb <- all_tox_wb %>%
  left_join(ENS_names[, c("gene_id", "gene_name")], by = "gene_id")
all_tox_wb$gene_name <- ifelse(is.na(all_tox_wb$gene_name.x), all_tox_wb$gene_name.y, all_tox_wb$gene_name.x)
df_wb_beforesens <- all_tox_wb %>% dplyr::select(-gene_name.x, -gene_name.y) #1403
str(df_wb_beforesens)
df_wb_beforesens$chromosome <- as.integer(df_wb_beforesens$chromosome)

wb_genes_tox_names_beforesens <- unique(df_wb_beforesens$gene_name) #313 (without NA)


####################################################################################################################################################################
####################load data from SENSITIVITY ANALYSIS########SCRIPT 4##############################################################################################
####################################################################################################################################################################

#go to different folder

random_all <- NULL
nrow_total <- 0

for (chemical in chemical_list){
  file_name <- paste0("sensitivity_analysis_4km/RANDOM_all_snps_wb_4km_",chemical,".txt")
  RANDOM_output <- read.table(file_name, header = T)
  n <- nrow(RANDOM_output)
  nrow_total <- nrow_total + n 
  RANDOM_output$tox_name <- chemical
  random_all <-rbind(random_all,RANDOM_output)
}

####################################################################################################################################################################
####################filter data from SENSITIVITY ANALYSIS############################################################################################################
####################################################################################################################################################################

#all rows
df_sa <- random_all #328,792 (100 runs per chemical)
df_sa <- random_all[c(1:3, 5)] #without qvalue

#only rows that occur min 2 times
df_sa_duplicates <- df_sa[duplicated(df_sa) | duplicated(df_sa, fromLast = TRUE), ]
df_sa_duplicates #191,458

#only rows that occur min 5 times
df_sa_min_5 <- df_sa %>%
  group_by(across(everything())) %>%  # Group by all columns
  filter(n() >= 5) %>%               # Filter groups with at least 5 occurrences
  ungroup()                          # Ungroup the data
#56,240

#only rows that occur min 10 times
df_sa_min_10 <- df_sa %>%
  group_by(across(everything())) %>%  # Group by all columns
  filter(n() >= 10) %>%               # Filter groups with at least 10 occurrences
  ungroup()  
#3081

###############################################################################################################################
#REMOVE SNPS FROM RDA output (script 3) based upon Sensitivity Analysis
###############################################################################################################################

#same variables, same order = snp, chromosome, position, tox_name
all_tox_wb_without_q <- df_wb_beforesens[c(1:3,8)] #4 variables (snp, chrom, pos, toxname) = 1403

#REMOVE ALL SNPS PICKED UP BY SA
df_after_sens <- dplyr::anti_join(all_tox_wb_without_q, df_sa) #827 rows 

#REMOVE SNPS THAT POP UP AT LEAST 2 TIMES
df_after_sens_duplicates <- dplyr::anti_join(all_tox_wb_without_q, df_sa_duplicates) #1238

#REMOVE SNPS THAT POP UP AT LEAST 5 TIMES
df_after_sens_min_5 <- dplyr::anti_join(all_tox_wb_without_q, df_sa_min_5) #1389
#REMOVE SNPS THAT POP UP AT LEAST 10 TIMES
df_after_sens_min_10 <- dplyr::anti_join(all_tox_wb_without_q, df_sa_min_10) #1403


###############################################################################################################################
#Use the filtered dataset with 4 variables and extend it to include gene names
###############################################################################################################################

df_after_sens_2 <- left_join(df_after_sens, df_wb_beforesens, by = c("snp","chromosome","position","tox_name"),relationship = "many-to-many")
df_after_sens_3 <- df_after_sens_2 %>% distinct()   #827 unique rows

removed_links <- df_wb_beforesens %>%
  anti_join(df_after_sens, by = c("snp","chromosome","position","tox_name")) #576

df_after_sens_duplicates_2 <- left_join(df_after_sens_duplicates, df_wb_beforesens, by = c("snp","chromosome","position","tox_name"),relationship = "many-to-many")
df_after_sens_duplicates_3 <- df_after_sens_duplicates_2 %>% distinct()   #1,238 unique rows

df_after_sens_min_5_2 <- left_join(df_after_sens_min_5, df_wb_beforesens, by = c("snp","chromosome","position","tox_name"),relationship = "many-to-many")
df_after_sens_min_5_3 <- df_after_sens_min_5_2 %>% distinct()    #1389

df_after_sens_min_10_2 <- left_join(df_after_sens_min_10, df_wb_beforesens, by = c("snp","chromosome","position","tox_name"),relationship = "many-to-many")
df_after_sens_min_10_3 <- df_after_sens_min_10_2 %>% distinct()    #1403

###############################################################################################################################
#Output
###############################################################################################################################

saveRDS(df_after_sens_3, "all_tox_wb_after_sens_4km.RData") 
saveRDS(removed_links, "removed_links_4km.RData") 
saveRDS(df_after_sens_duplicates_3, "all_tox_wb_after_sens_duplicate_4km.RData") 
saveRDS(df_after_sens_min_5_3, "all_tox_wb_after_sens_min_5_4km.RData") 
saveRDS(df_after_sens_min_10_3, "all_tox_wb_after_sens_min_10_4km.RData") 



unique(df_after_sens_3$tox_name) #49
filtered_df <- df_after_sens_3[is.na(df_after_sens_3$gene_name), ] #480
filtered_df_2 <- df_after_sens_3[!is.na(df_after_sens_3$gene_name), ] #347
unique(filtered_df_2$tox_name) #38
unique(filtered_df_2$gene_name) #212

tox_gene <- unique(filtered_df_2[c(4,9)]) #325
df_count_chem <- tox_gene %>%group_by(tox_name) %>%summarise(count = n()) #38
df_count_gene <- tox_gene %>%group_by(gene_name) %>%summarise(count = n()) #212
