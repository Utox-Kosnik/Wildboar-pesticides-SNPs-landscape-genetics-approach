############################################################################################
###################### CHECKING EFFECT OF COUNTRY #########################################
############################################################################################

Sys.setenv(LANG = "en")
rm(list = ls()) 

######################
#Packages needed 
######################

library(dplyr)
library(data.table)
library(ggplot2)

######################
#Input files needed
######################

#MATRIX FILTERED FOR 282 INDIVIDUALS
wb_matrix_no_na_282 <- readRDS("wb_matrix_no_na_282.RData") #FROM SCRIPT 2

#FINAL METADATA
wb_final_metadata <-readRDS("wb_final_metadata_4km.RData") #FROM SCRIPT 2

#OUTPUT
df_after_sens <- readRDS("all_tox_wb_after_sens_duplicate_4km.RData") #FROM SCRIPT 5
wb_genes_tox <- df_after_sens %>% dplyr::filter(!is.na(gene_name)) 
unique(wb_genes_tox$gene_name) #212 (singletons removed) vs 278 (duplicates removed)

top5 <- as.data.frame(sort(table(wb_genes_tox$gene_name),decreasing = TRUE)[1:5])
top5
#PDE4D (9)
#RASSF5 (9)
#KCNMA1 (8)
#PLCXD3 (8)
#GPC5 (7)

#pick a gene
df_1  <- wb_genes_tox %>% filter(gene_name == "PDE4D")
df_2 <- wb_genes_tox %>% filter(gene_name == "RASSF5")
df_3 <- wb_genes_tox %>% filter(gene_name == "KCNMA1")
df_4 <- wb_genes_tox %>% filter(gene_name == "PLCXD3")
df_5 <- wb_genes_tox %>% filter(gene_name == "GPC5")

wb_id <- rownames(wb_matrix_no_na_282)
country <- as.data.frame(wb_final_metadata[,c(6)])

####################################################################################################################################################################
############CHECK INDIVIDUAL GENES########################################################################################################################################################
####################################################################################################################################################################

#PDE4D
#e.g.,: SNP names INRA0051559, MARC0084747, DRGA0016146, ALGA0090426, ASGA0073164 
#SNP nrs  37605, 37608,  37611, 37615, 37618 (check names_SNP.txt file)

m_pde4d <- wb_matrix_no_na_282[,c(37615)] 
setDT(df_1) #convert to data table
filtered_df_pde4d <- df_1[snp == "ALGA0090426"]
df_summary = NULL

for (chemical in filtered_df_pde4d$tox_name){
toxic_column <- unlist(wb_final_metadata[chemical]) 
new_df <- as.data.frame(cbind(wb_id, m_pde4d, toxic_column, country))
colnames(new_df) <- c("id", "m", "tox", "country")
new_df$tox <- as.numeric(new_df$tox)
new_df$m <- factor(new_df$m,levels = c("0", "1", "2"),ordered = TRUE)

new_df <- new_df[!is.na(new_df$tox), ]
country_list <- unique(new_df$country)
for (c in country_list){
  country_data <- subset(new_df, country == c)
  if (length(unique(country_data$m)) > 1){
    kr <- kruskal.test(tox ~ m, data = country_data)
    if (!is.na(kr$p.value) && kr$p.value < 0.05){
      chem_country_combo <- cbind(chemical,c)
      df_summary <- rbind(df_summary, chem_country_combo)
      }
  }
}}

ggplot(data=new_df, aes(x=m, y=tox, color= country)) + geom_boxplot() + theme_minimal()#

####################################################################################################################################################################

#RASSF5
#SNP name	ALGA0056121
#SNP nr 23033

m_rassf5 <- wb_matrix_no_na_282[,c(23033)] 
setDT(df_2) #convert to data table
filtered_df_rassf5 <- df_2[snp == "ALGA0056121"]
df_summary = NULL

for (chemical in filtered_df_rassf5$tox_name){
toxic_column <- unlist(wb_final_metadata[chemical]) 
new_df <- as.data.frame(cbind(wb_id, m_rassf5, toxic_column, country))
colnames(new_df) <- c("id", "m", "tox", "country")
new_df$tox <- as.numeric(new_df$tox)
new_df$m <- factor(new_df$m,levels = c("0", "1", "2"),ordered = TRUE)

new_df <- new_df[!is.na(new_df$tox), ]
country_list <- unique(new_df$country)
for (c in country_list){
  country_data <- subset(new_df, country == c)
  if (length(unique(country_data$m)) > 1){
    kr <- kruskal.test(tox ~ m, data = country_data)
    if (!is.na(kr$p.value) && kr$p.value < 0.05){
      chem_country_combo <- cbind(chemical,c)
      df_summary <- rbind(df_summary, chem_country_combo)}
  }
}}

ggplot(data=new_df, aes(x=m, y=tox, color= country)) + geom_boxplot()  

####################################################################################################################################################################

#KCNMA1
#e.g., SNP name INRA0045133, ALGA0079181
#SNP nr 33239, 33243

m_kcnma1 <- wb_matrix_no_na_282[,c(33239)] 
setDT(df_3) #convert to data table
filtered_df_kcnma1 <- df_3[snp == "INRA0045133"]
df_summary = NULL

for (chemical in filtered_df_kcnma1$tox_name){
toxic_column <- unlist(wb_final_metadata[chemical]) 
new_df <- as.data.frame(cbind(wb_id, m_kcnma1, toxic_column, country))
colnames(new_df) <- c("id", "m", "tox", "country")
new_df$tox <- as.numeric(new_df$tox)
new_df$m <- factor(new_df$m,levels = c("0", "1", "2"),ordered = TRUE)

new_df <- new_df[!is.na(new_df$tox), ]
country_list <- unique(new_df$country)
for (c in country_list){
  country_data <- subset(new_df, country == c)
  if (length(unique(country_data$m)) > 1){
    kr <- kruskal.test(tox ~ m, data = country_data)
    if (!is.na(kr$p.value) && kr$p.value < 0.05){
      chem_country_combo <- cbind(chemical,c)
      df_summary <- rbind(df_summary, chem_country_combo)
    }
  }
}}

ggplot(data=new_df, aes(x=m, y=tox, color= country)) + geom_boxplot() + geom_jitter() 

####################################################################################################################################################################

#PLCXD3
#e.g., SNP name ALGA0089898
#SNP nr 37347

m_plcxd3 <- wb_matrix_no_na_282[,c(37347)] 
setDT(df_4) #convert to data table
filtered_df_plcxd3 <- df_4[snp == "ALGA0089898"]
df_summary = NULL

for (chemical in filtered_df_plcxd3$tox_name){
toxic_column <- unlist(wb_final_metadata[chemical]) 
new_df <- as.data.frame(cbind(wb_id, m_plcxd3, toxic_column, country))
colnames(new_df) <- c("id", "m", "tox", "country")
new_df$tox <- as.numeric(new_df$tox)
new_df$m <- factor(new_df$m,levels = c("0", "1", "2"),ordered = TRUE)
new_df <- new_df[!is.na(new_df$tox), ]
country_list <- unique(new_df$country)
for (c in country_list){
  country_data <- subset(new_df, country == c)
  if (length(unique(country_data$m)) > 1){
    kr <- kruskal.test(tox ~ m, data = country_data)
    if (!is.na(kr$p.value) && kr$p.value < 0.05){
      chem_country_combo <- cbind(chemical,c)
      df_summary <- rbind(df_summary, chem_country_combo)
    }
  }
}}

ggplot(data=new_df, aes(x=m, y=tox, color= country)) + geom_boxplot() + geom_jitter() 

####################################################################################################################################################################

#GPC5
#e.g., SNP names DRGA0011366, MARC0054697
#SNP nr 26580, 26579

m_gpc5 <- wb_matrix_no_na_282[,c(26580)] 
setDT(df_5) #convert to data table
filtered_df_gpc5 <- df_5[snp == "DRGA0011366"]
df_summary = NULL

for (chemical in filtered_df_gpc5$tox_name){
print(chemical)
toxic_column <- unlist(wb_final_metadata[chemical]) 
new_df <- as.data.frame(cbind(wb_id, m_gpc5, toxic_column, country))
colnames(new_df) <- c("id", "m", "tox", "country")
new_df$tox <- as.numeric(new_df$tox)
new_df$m <- factor(new_df$m,levels = c("0", "1", "2"),ordered = TRUE)

new_df <- new_df[!is.na(new_df$tox), ]
country_list <- unique(new_df$country)
for (c in country_list){
  country_data <- subset(new_df, country == c)
  if (length(unique(country_data$m)) > 1){
    kr <- kruskal.test(tox ~ m, data = country_data)
    if (!is.na(kr$p.value) && kr$p.value < 0.05){
      print(c)
      chem_country_combo <- cbind(chemical,c)
      df_summary <- rbind(df_summary, chem_country_combo)
    }
  }
}}

ggplot(data=new_df, aes(x=m, y=tox, color= country)) + geom_boxplot() + geom_jitter() 
