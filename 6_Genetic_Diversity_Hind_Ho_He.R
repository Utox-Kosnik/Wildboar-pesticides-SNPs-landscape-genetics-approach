############################################################################################
###################### GENETIC DIVERSITY ANALYSIS #########################################
############################################################################################

Sys.setenv(LANG = "en")
rm(list = ls()) 

######################
#Packages needed (some packages are called within the code)
######################

library(dplyr) #%>%
library(rsq)
library(ggplot2)


test_match_order2 <- function(x,y) {
  if (isTRUE(all.equal(x,y))) print('Perfect match in same order')
  if (!isTRUE(all.equal(x,y)) && isTRUE(all.equal(sort(x),sort(y)))) print('Perfect match in wrong order')
  if (!isTRUE(all.equal(x,y)) && !isTRUE(all.equal(sort(x),sort(y)))) print('No match')
}

######################
#Input files
######################

#CHEMICAL_LIST
chemical_list <- readRDS("chemical_list.RData") #FROM SCRIPT 2

#MATRIX FILTERED FOR 282 INDIVIDUALS
wb_matrix_no_na_282 <- readRDS("wb_matrix_no_na_282.RData") #FROM SCRIPT 2

#FINAL METADATA
wb_final_metadata <-readRDS("wb_final_metadata_4km.RData") #FROM SCRIPT 2

#METADATA FOR POPULATION LEVEL ANALYSES
df_pop <-  read.table("input_file_for_population_analyses.txt", header = T) #only individuals with min 4 ind per same location
df_pop$pop_nr <- as.factor(df_pop$location)
df_pop <- df_pop %>% filter(country != "Serbia" & country != "Bosnia") #BECAUSE THESE DO NOT HAVE pesticide data

df_pop_final <- merge(df_pop,wb_final_metadata,by="id")
wb_matrix_no_na_199 <- wb_matrix_no_na_282[rownames(wb_matrix_no_na_282) %in% df_pop_final$id, ] #keep 199 individuals
test_match_order2(rownames(wb_matrix_no_na_199), df_pop_final$id) # "Perfect match in wrong order"
right_order <- rownames(wb_matrix_no_na_199)
df_pop_final <- df_pop_final[match(right_order, df_pop_final$id),]
test_match_order2(rownames(wb_matrix_no_na_199), df_pop_final$id)

final_pest_per_pop <- df_pop_final %>% group_by(pop_nr) %>% dplyr::slice(1) %>% ungroup()

######################################################################################################################################
#INDIVIDUAL-LEVEL HETEROZYGOSITY
######################################################################################################################################

number_heterozyg_loci <- rowSums(wb_matrix_no_na_282 == 1)
heterozygosity = (number_heterozyg_loci / 40700)
df_hzyg <- cbind(heterozygosity, wb_final_metadata)

min(df_hzyg$heterozygosity) #10% = W-Sardinia
max(df_hzyg$heterozygosity) #45% = Estonia
mean(df_hzyg$heterozygosity) #23%
median(df_hzyg$heterozygosity) #23%

kruskal.test(heterozygosity ~ country, data = df_hzyg) #p<0.05
ggplot(df_hzyg, aes(x = country, y = heterozygosity)) + geom_boxplot()

i = 0
significant_chemicals <- c()
for (chemical in chemical_list){
  pest <- df_hzyg[[chemical]]
  df_hzyg_no_na <- df_hzyg %>% filter(!is.na(pest))
  df_hzyg_no_na$country <- droplevels(df_hzyg_no_na$country)
  pest_no_na <- na.omit(pest)
  
  if(nlevels(df_hzyg_no_na$country) > 1) {

    model1 <- glm(heterozygosity ~ pest_no_na + latitude + longitude + agricultural_percentage + human_population, family=quasibinomial, data = df_hzyg_no_na)
    model2 <- glm(heterozygosity ~ latitude  + longitude + agricultural_percentage + human_population, family=quasibinomial, data = df_hzyg_no_na)
    an <- anova(model1, model2, test="Chisq")
    quasi_r2_model1 <- rsq(model1)
    quasi_r2_model2 <- rsq(model2)
    R2_difference <- quasi_r2_model1 - quasi_r2_model2
    
  }
  if (an$`Pr(>Chi)`[2] < 0.05) { #if pvalue <0.05
    
    estimate <- model1$coefficients[2]
    result <- cbind(chemical, estimate, R2_difference)
    print(result)
    i = i + 1 #40
    significant_chemicals <- rbind(significant_chemicals, result)
  }
}

#40 chemicals significantly explain variation in heterozygosity

hzyg_results <- as.data.frame(significant_chemicals)
rownames(hzyg_results) <- c(1:i)
hzyg_results$estimate <- as.numeric(hzyg_results$estimate)
hzyg_results$R2_difference <- as.numeric(hzyg_results$R2_difference)

min(hzyg_results$R2_difference) #0.09%
max(hzyg_results$R2_difference) #21%
mean(hzyg_results$R2_difference) #3.5%
median(hzyg_results$R2_difference) #2.2%

write.table(hzyg_results, "heterozygosity_individual.txt")

################################################################################################################################
#Observed/Expected heterozygosity 
###############################################################################################################################

#Observed and expected heterozygosity
population_list <- as.numeric(unique(df_pop_final$pop_nr))

obs_heteroz_all_pop = NULL
exp_heteroz_all_pop = NULL
freq_all_1 = 0
freq_all_2 = 0

for (p in population_list){ #loop through populations
  print(p)
  df_pop_nr <- df_pop_final %>% filter(pop_nr == p)
  ids_to_keep_pop <- df_pop_nr$id
  n_total <- length(ids_to_keep_pop)
  wb_matrix_no_na_pop <- wb_matrix_no_na_199[rownames(wb_matrix_no_na_199) %in% ids_to_keep_pop, ] 
  
  #observed heterozygosity
  n_heteroz_counts <- colSums(wb_matrix_no_na_pop == 1) #number of heterozygotes per SNP
  obs_heteroz_per_snp <- n_heteroz_counts/n_total #obs heteroz percentage per SNP
  obs_heteroz_1_pop <- sum(obs_heteroz_per_snp)/40700 #obs heteroz percentage for all SNPs
  
  #expected heterozygosity
  n_homoz_counts_all_1 <- colSums(wb_matrix_no_na_pop == 0) #number of homozygotes for allele 1 per SNP
  n_homoz_counts_all_2 <- colSums(wb_matrix_no_na_pop == 2) #number of homozygotes for allele 2 per SNP
  freq_all_1 <- ((n_homoz_counts_all_1*2)+n_heteroz_counts)/(n_total*2) #freq allele 1 per SNP
  freq_all_2 <- ((n_homoz_counts_all_2*2)+n_heteroz_counts)/(n_total*2)  #freq allele 2 per SNP
  exp_heteroz_per_snp <- 1-(freq_all_1^2)-(freq_all_2^2) #exp heterozyg per SNP
  exp_heteroz_1_pop <- sum(exp_heteroz_per_snp)/40700 #exp heterozyg for all SNPs
  
  #
  obs_heteroz_all_pop <- cbind(obs_heteroz_all_pop, obs_heteroz_1_pop)
  exp_heteroz_all_pop <- cbind(exp_heteroz_all_pop, exp_heteroz_1_pop)
}
colnames(obs_heteroz_all_pop) <- population_list
colnames(exp_heteroz_all_pop) <- population_list
heteroz_all_pop <- rbind(obs_heteroz_all_pop, exp_heteroz_all_pop)
rownames(heteroz_all_pop)<-c("Ho","He")

heteroz_all_pop_t <- as.data.frame(t(heteroz_all_pop))
heteroz_all_pop_t$pop_nr <- rownames(heteroz_all_pop_t)
mean(heteroz_all_pop_t$Ho) #23%
mean(heteroz_all_pop_t$He) #21%
min(heteroz_all_pop_t$Ho) #15%, West-Sardinia
min(heteroz_all_pop_t$He) #14%, West-Sardinia
max(heteroz_all_pop_t$Ho) #35%, Estonia
max(heteroz_all_pop_t$He) #35%, Estonia

final_HO_HE <- merge(final_pest_per_pop, heteroz_all_pop_t,by="pop_nr") #every population has values for Ho, He and pesticide concentrations

#Observed H

i = 0
check = 0
significant_chemicals_Ho <- c()
for (chemical in chemical_list){
  pest <- final_HO_HE[[chemical]]
  final_HO_HE_no_na <- final_HO_HE %>% filter(!is.na(pest))
  pest_no_na <- na.omit(pest)
  if (nrow(final_HO_HE_no_na) > 5 ){ #min 5 populations to test
    check = check + 1 # chemicals checked = 116
    model1 <- lm(log(Ho) ~ pest_no_na + latitude.x + longitude.x + agricultural_percentage + human_population, data = final_HO_HE_no_na)
    model2 <- lm(log(Ho) ~ latitude.x  + longitude.x + agricultural_percentage + human_population, data = final_HO_HE_no_na)
    an <- anova(model1, model2, test="Chisq")
    summary(model1)
    m1_adjusted_r_squared <- summary(model1)$adj.r.squared
    m2_adjusted_r_squared <- summary(model2)$adj.r.squared
    
    if (m1_adjusted_r_squared < 0){
      m1_adjusted_r_squared <- 0
    }
    
    if (m2_adjusted_r_squared < 0){
      m2_adjusted_r_squared <- 0
    }
    
    R2_difference <- m1_adjusted_r_squared - m2_adjusted_r_squared
    
    if (an$`Pr(>Chi)`[2] < 0.05) { #if the p-value of the anova is lower than 0.05
      estimate <- model1$coefficients[2]
      result <- cbind(chemical, estimate, R2_difference, nrow(final_HO_HE_no_na))
      print(result)
      i = i + 1 #
      # If significant, save the column name
      significant_chemicals_Ho <- rbind(significant_chemicals_Ho, result)
    }}
}

#Ho : 5 chemicals significantly explain variation in observed heterozygosity 
colnames(significant_chemicals_Ho)[4] <- c("populations_included")
rownames(significant_chemicals_Ho) <- c(1:i)

write.table(significant_chemicals_Ho, "heterozygosity_observed.txt")

#Expected H

i = 0
check = 0
significant_chemicals_He <- c()
for (chemical in chemical_list){
  pest <- final_HO_HE[[chemical]]
  final_HO_HE_no_na <- final_HO_HE %>% filter(!is.na(pest))
  pest_no_na <- na.omit(pest)
  if (nrow(final_HO_HE_no_na) > 5 ){ #min 5 populations to test
    check = check + 1 # chemicals checked = 116
    model1 <- lm(log(He) ~ pest_no_na + latitude.x + longitude.x + agricultural_percentage + human_population, data = final_HO_HE_no_na)
    model2 <- lm(log(He) ~ latitude.x  + longitude.x + agricultural_percentage + human_population, data = final_HO_HE_no_na)
    an <- anova(model1, model2, test="Chisq")
    summary(model1)
    m1_adjusted_r_squared <- summary(model1)$adj.r.squared
    m2_adjusted_r_squared <- summary(model2)$adj.r.squared
    
    if (m1_adjusted_r_squared < 0){
      m1_adjusted_r_squared <- 0
    }
    
    if (m2_adjusted_r_squared < 0){
      m2_adjusted_r_squared <- 0
    }
    
    R2_difference <- m1_adjusted_r_squared - m2_adjusted_r_squared
    
    if (an$`Pr(>Chi)`[2] < 0.05) { #if the p-value of the anova is lower than 0.05
      estimate <- model1$coefficients[2]
      result <- cbind(chemical, estimate, R2_difference, nrow(final_HO_HE_no_na))
      print(result)
      i = i + 1 #
      # If significant, save the column name
      significant_chemicals_He <- rbind(significant_chemicals_He, result)
    }}
}

#He : 7 chemicals significantly explain variation in expected heterozygosity (positive relationship)
colnames(significant_chemicals_He)[4] <- c("populations_included")
rownames(significant_chemicals_He) <- c(1:i)

write.table(significant_chemicals_He, "heterozygosity_expected.txt")
