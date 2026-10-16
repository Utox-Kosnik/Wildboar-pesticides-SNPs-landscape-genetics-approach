#################################################################################################################################
################################# IDENTIFY ALL GENES PRESENT IN THE INPUT DATA (40,700 SNPs) #####################################
#################################################################################################################################

Sys.setenv(LANG = "en")
rm(list = ls()) 

library(dplyr) #%>%
######################
#Input files needed
######################

#SNP names, taken from map file
all_snp_names_from_map <- read.table("names_SNPs.txt", header = F) #FROM MAP file

#ILLUMINA DOCUMENT:match SNP name with real chromosome + position (to later version of the reference genome: from 10.2 to 11.1)
snp_id<- read.table("forR_PorcineSNP60.csv",  sep=",", header = T) #FROM ILLUMINA

#GENOME DATA TO MATCH SNPS WITH GENES
g <- rtracklayer::readGFF("Ensembl_Sus_scrofa.Sscrofa11.1.112.gtf") #FROM ENSEMBL
str(g)
g$start<-as.numeric(g$start)
g$end <- as.numeric(g$end)

#######################################################
#1: transfer from reference genome 10.2 to 11.1
#########################################################

colnames(all_snp_names_from_map)<-c("chrom","name","x","pos")
all_snp_test<- all_snp_names_from_map
all_snp_test<- merge(all_snp_test, snp_id[, c("Name", "Chr", "MapInfo")], by.x="name", by.y="Name")
all_snp_test<- data.table::as.data.table(all_snp_test)
all_snp <- all_snp_test[, list(snp=name, chromosome=Chr, position=MapInfo)]
str(all_snp)
all_snp$position <- as.numeric(all_snp$position)

#######################################################################################################
#MATCH SNPS TO GENES
####################################################################################################### 

final_df = NULL
n=0
for(i in 1:nrow(all_snp)){
  print(i)
  right_snp <- all_snp$snp[i]
  right_chromosome <- all_snp$chr[i]
  right_position <- all_snp$position[i]
  g_filter <- g %>% dplyr::filter(seqid == right_chromosome, type == "gene", 
                                    (start <= right_position & right_position <= end)) %>% as.data.frame()  #intragenic SNPs
  g_filter <- g_filter[!duplicated(g_filter), ]
  if (nrow(g_filter) == 0){ #no closeby genes
    g_subset <- g_filter[nrow(g_filter)+1,] <- NA
    g_subset <- g_filter[,c(1:5,9,11)]
  } else { #intragenic genes
    n <- n + 1 #to know how many snps are located within a gene
    g_subset <-g_filter[,c(1:5,9,11)] #index by name
  }
  g_subset$snp <- right_snp
  g_subset$chromosome <- right_chromosome
  g_subset$position <- right_position
  final_df = rbind(final_df,g_subset)
  }

saveRDS(final_df, "all_genes_represented_in_40700_inputdata.RData") 
sum(!is.na(final_df$gene_name)) #17,240 SNPs are intragenic --> 17,643 genes
sum(is.na(final_df$gene_name)) #23,460 SNPs are intergenic
length(unique(na.omit(final_df$gene_name))) #6551 unique genes

