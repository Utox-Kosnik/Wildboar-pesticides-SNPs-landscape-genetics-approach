############################################################################
#fix missing genotypes WILD BOAR
###########################################################################

rm(list = ls()) 

#Packages needed 

library(adegenet)
library(dplyr)

packageVersion("adegenet")
citation()

#Input files needed

wb <- read.PLINK("wb_eur_0.1_maf_0.01.raw")
wb_matrix <- as.matrix(wb, "matrix") #306 individuals, 40700 SNPs
sum(is.na(wb_matrix)) #10,839 NAs in dataset
metadata_simple_wb <- read.table("metadata_simple_wb.txt", header = T) 

test_match_order2 <- function(x,y) {
  
  if (isTRUE(all.equal(x,y))) print('Perfect match in same order')
  
  if (!isTRUE(all.equal(x,y)) && isTRUE(all.equal(sort(x),sort(y)))) print('Perfect match in wrong order')
  
  if (!isTRUE(all.equal(x,y)) && !isTRUE(all.equal(sort(x),sort(y)))) print('No match')
}

###########################################################################

#check whether the rownames of the final metadata and final matrix are the same
test_match_order2(rownames(wb_matrix), (metadata_simple_wb$id)) #"Perfect match in wrong order"
right_order <- rownames(wb_matrix)
metadata_simple_wb <- metadata_simple_wb[match(right_order, metadata_simple_wb$id),]
test_match_order2(rownames(wb_matrix), (metadata_simple_wb$id)) #"Perfect match in same order"

# Create a new matrix to store imputed values
wb_matrix_noNA <- wb_matrix

# Find positions of NAs
na_positions <- which(is.na(wb_matrix), arr.ind = TRUE)

for (pos in 1:nrow(na_positions)) {
  row_idx <- na_positions[pos, 1]  # Individual (row)
  col_idx <- na_positions[pos, 2]  # SNP (column)
  
  individual_id <- rownames(wb_matrix)[row_idx]
  
  # Find corresponding pca_label
  pop_name <- metadata_simple_wb$pop_name[metadata_simple_wb$id == individual_id]
  
  # Get all individuals with the same pca_label
  same_group <- metadata_simple_wb$id[metadata_simple_wb$pop_name == pop_name]
  
  # Extract the values for the same SNP (column) within the group
  group_values <- wb_matrix_noNA[same_group, col_idx, drop = FALSE]
  
  # Remove NAs from the group
  group_values <- group_values[!is.na(group_values)]
  
  # Check if the group has valid values
  if (length(group_values) > 0) {
    # Find the most frequent allele within the group
    allele_counts <- table(group_values)
    most_frequent_values <- as.numeric(names(allele_counts)[allele_counts == max(allele_counts)])
    
    # If there's a tie, assign 1
    if (length(most_frequent_values) > 1) {
      most_frequent <- 1
    } else {
      most_frequent <- most_frequent_values
    }
    
  } else {
    # If group is empty (all NAs), check all individuals in the dataset
    all_values <- wb_matrix_noNA[, col_idx]
    all_values <- all_values[!is.na(all_values)]  # Remove NAs
    
    # Only proceed if there are valid values in the dataset
    if (length(all_values) > 0) {
      allele_counts <- table(all_values)
      most_frequent_values <- as.numeric(names(allele_counts)[allele_counts == max(allele_counts)])
      
      # If there's a tie, assign 1
      if (length(most_frequent_values) > 1) {
        most_frequent <- 1
      } else {
        most_frequent <- most_frequent_values
      }
    } else {
      # If the entire dataset has NAs for this SNP, keep NA
      next
    }
  }
  
  # Assign the imputed value in the new matrix
  wb_matrix_noNA[row_idx, col_idx] <- most_frequent
}
storage.mode(wb_matrix_noNA) <- "integer"
non_na_indices <- !is.na(wb_matrix)
identical(wb_matrix[non_na_indices], wb_matrix_noNA[non_na_indices]) #TRUE

saveRDS(wb_matrix_noNA, "final_matrix_wb_europe_no_na.RData")
