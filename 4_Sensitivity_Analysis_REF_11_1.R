############################################################################################
################################# SENSITIVITY ANALYSIS #####################################
############################################################################################

Sys.setenv(LANG = "en")
rm(list = ls())


############################################################################################
# PACKAGES
############################################################################################

library(dplyr)

# Other packages are called explicitly in the code:
# vegan
# qvalue
# tidyr


############################################################################################
# INPUT FILES
############################################################################################

# CHEMICAL LIST
chemical_list <- readRDS("chemical_list.RData")   # FROM SCRIPT 2


# MATRIX FILTERED FOR 282 INDIVIDUALS
wb_matrix_no_na_282 <- readRDS(
  "wb_matrix_no_na_282.RData"
)   # FROM SCRIPT 2


# FINAL METADATA - 4 KM
wb_final_metadata <- readRDS(
  "wb_final_metadata_4km.RData"
)   # FROM SCRIPT 2


# ILLUMINA DOCUMENT:
# Match SNP name with chromosome + position in Sus scrofa reference genome 11.1
snp_id <- read.table(
  "forR_PorcineSNP60.csv",
  sep = ",",
  header = TRUE
)


# SNP names taken from MAP file
all_snp_names_from_map <- read.table(
  "names_SNPs.txt",
  header = FALSE
)


############################################################################################
# CHECK SAMPLE ORDER
############################################################################################

test_match_order2 <- function(x, y) {
  
  if (isTRUE(all.equal(x, y))) {
    print("Perfect match in same order")
  }
  
  if (
    !isTRUE(all.equal(x, y)) &&
    isTRUE(all.equal(sort(x), sort(y)))
  ) {
    print("Perfect match in wrong order")
  }
  
  if (
    !isTRUE(all.equal(x, y)) &&
    !isTRUE(all.equal(sort(x), sort(y)))
  ) {
    print("No match")
  }
}


test_match_order2(
  rownames(wb_matrix_no_na_282),
  wb_final_metadata$id
)


############################################################################################
# PREPARE MAP FILE
############################################################################################

# The MAP file chromosome and position columns refer to an older reference genome.
# Therefore, these columns are only used here to identify the SNP names.
# The correct chromosome and position are later obtained from snp_id.

all_snp_names_from_map$nr <- rownames(all_snp_names_from_map)

colnames(all_snp_names_from_map) <- c(
  "chrom",
  "name",
  "x",
  "pos",
  "nr"
)


############################################################################################
############################################################################################
############################## SENSITIVITY ANALYSIS ########################################
############################################################################################
############################################################################################


for (chemical in chemical_list) {
  
  
  ##########################################################################################
  # INITIALISE OUTPUT FOR THIS CHEMICAL
  ##########################################################################################
  
  all_snps_wb <- NULL
  
  failed_seeds <- integer(0)
  
  
  print("------------------------------------------------------------")
  print(paste("Chemical:", chemical))
  
  
  ##########################################################################################
  # EXPOSURE COLUMN
  ##########################################################################################
  
  toxic_column <- unlist(
    wb_final_metadata[chemical]
  )
  
  
  print(
    paste(
      "Number of unique exposure values:",
      length(unique(toxic_column))
    )
  )
  
  
  ##########################################################################################
  # DEFINE NA POSITIONS ONCE
  ##########################################################################################
  
  na_positions <- is.na(toxic_column)
  
  non_na_values <- toxic_column[!na_positions]
  
  
  ##########################################################################################
  # RUN 100 SUCCESSFUL RANDOM PERMUTATIONS
  ##########################################################################################
  
  target_success <- 100
  
  n_success <- 0
  
  jjj <- 1
  
  
  while (n_success < target_success) {
    
    
    print(
      paste0(
        "Chemical = ", chemical,
        " | seed = ", jjj,
        " | successful permutation = ", n_success + 1,
        "/", target_success
      )
    )
    
    
    ########################################################################################
    # RANDOMISE EXPOSURE VALUES
    ########################################################################################
    
    # Always start from the ORIGINAL exposure column
    shuffled_toxic_column <- toxic_column
    
    
    set.seed(jjj)
    
    
    shuffled_toxic_column[!na_positions] <- sample(
      non_na_values
    )
    
    
    ########################################################################################
    # RDA
    #
    # Occasionally LAPACK's SVD algorithm may fail for a particular random permutation.
    # tryCatch() prevents the entire sensitivity analysis from stopping.
    ########################################################################################
    
    wb.rda.tox <- tryCatch(
      
      {
        
        vegan::rda(
          
          wb_matrix_no_na_282 ~
            shuffled_toxic_column +
            Condition(
              longitude +
                latitude +
                agricultural_percentage +
                human_population
            ),
          
          data = wb_final_metadata,
          
          scale = TRUE,
          
          na.action = na.omit
          
        )
        
      },
      
      error = function(e) {
        
        message(
          "RDA FAILED for ",
          chemical,
          " at seed ",
          jjj,
          ": ",
          conditionMessage(e)
        )
        
        return(NULL)
        
      }
      
    )
    
    
    ########################################################################################
    # IF RDA FAILED: SKIP THIS SEED
    ########################################################################################
    
    if (is.null(wb.rda.tox)) {
      
      failed_seeds <- c(
        failed_seeds,
        jjj
      )
      
      jjj <- jjj + 1
      
      next
      
    }
    
    
    ########################################################################################
    # RDA LOADINGS
    ########################################################################################
    
    load.rda.class <- vegan::scores(
      wb.rda.tox,
      choices = 1,
      display = "species"
    )
    
    
    load.rda.class <- as.data.frame(
      load.rda.class
    )
    
    
    load.rda.class$locus <- rownames(
      load.rda.class
    )
    
    
    ########################################################################################
    # CONVERT LOADINGS TO LONG FORMAT
    ########################################################################################
    
    mload.rda.class <- tidyr::gather(
      load.rda.class,
      axis,
      loading,
      -locus
    )
    
    
    ########################################################################################
    # STANDARDISE LOADINGS
    ########################################################################################
    
    mload.rda.class$scaled_loading <- as.numeric(
      scale(
        mload.rda.class$loading
      )
    )
    
    
    ########################################################################################
    # GENOMIC INFLATION FACTOR
    ########################################################################################
    
    lambda.class <- median(
      mload.rda.class$scaled_loading^2,
      na.rm = TRUE
    ) / 0.456
    
    
    ########################################################################################
    # P-VALUES
    ########################################################################################
    
    mload.rda.class$pvalues <- (
      1 -
        pnorm(
          abs(
            mload.rda.class$scaled_loading
          )
        )
    ) * 2
    
    
    ########################################################################################
    # Q-VALUES
    ########################################################################################
    
    qvalues.class <- qvalue::qvalue(
      mload.rda.class$pvalues
    )
    
    
    mload.rda.class$qvalues <- qvalues.class$qvalues
    
    
    df_rda.qvalues.class <- as.data.frame(
      mload.rda.class
    )
    
    
    ########################################################################################
    # SIGNIFICANT SNPS
    ########################################################################################
    
    snp_wb_tox <- subset(
      df_rda.qvalues.class,
      qvalues < 0.05
    )
    
    
    ########################################################################################
    # ONLY CONTINUE SNP MATCHING IF SIGNIFICANT SNPS WERE FOUND
    ########################################################################################
    
    if (nrow(snp_wb_tox) > 0) {
      
      
      df_rda_final <- snp_wb_tox[
        ,
        c(
          "locus",
          "qvalues"
        )
      ]
      
      
      colnames(df_rda_final) <- c(
        "nr",
        "qvalue"
      )
      
      
      # Ensure matching format
      df_rda_final$nr <- as.character(
        df_rda_final$nr
      )
      
      
      ######################################################################################
      # MATCH LOCUS NUMBER TO SNP NAME FROM MAP FILE
      ######################################################################################
      
      snp_matches <- merge(
        
        df_rda_final,
        
        all_snp_names_from_map[
          ,
          c(
            "nr",
            "name"
          )
        ],
        
        by = "nr",
        
        all.x = TRUE,
        
        sort = FALSE
        
      )
      
      
      ######################################################################################
      # MATCH SNP NAME TO CHROMOSOME + POSITION IN REFERENCE GENOME 11.1
      ######################################################################################
      
      snp_matches <- merge(
        
        snp_matches,
        
        snp_id[
          ,
          c(
            "Name",
            "Chr",
            "MapInfo"
          )
        ],
        
        by.x = "name",
        
        by.y = "Name",
        
        all.x = TRUE,
        
        sort = FALSE
        
      )
      
      
      ######################################################################################
      # CREATE FINAL SNP TABLE FOR THIS PERMUTATION
      ######################################################################################
      
      all_snp <- data.frame(
        
        snp = snp_matches$name,
        
        chromosome = snp_matches$Chr,
        
        position = as.numeric(
          snp_matches$MapInfo
        ),
        
        qvalue = as.numeric(
          snp_matches$qvalue
        )
        
      )
      
      
      ######################################################################################
      # ADD TO RESULTS FOR THIS CHEMICAL
      ######################################################################################
      
      all_snps_wb <- rbind(
        all_snps_wb,
        all_snp
      )
      
    }
    
    
    ########################################################################################
    # PERMUTATION COMPLETED SUCCESSFULLY
    ########################################################################################
    
    n_success <- n_success + 1
    
    
    ########################################################################################
    # MOVE TO NEXT SEED
    ########################################################################################
    
    jjj <- jjj + 1
    
    
  }   # END WHILE LOOP
  
  
  ##########################################################################################
  # REPORT FAILED SEEDS
  ##########################################################################################
  
  print("------------------------------------------------------------")
  
  
  if (length(failed_seeds) > 0) {
    
    message(
      chemical,
      ": failed RDA seeds = ",
      paste(
        failed_seeds,
        collapse = ", "
      )
    )
    
  } else {
    
    message(
      chemical,
      ": no RDA permutations failed."
    )
    
  }
  
  
  print(
    paste(
      "Successful permutations:",
      n_success
    )
  )
  
  
  ##########################################################################################
  # IF NO SIGNIFICANT SNP WAS FOUND ACROSS ANY OF THE 100 PERMUTATIONS
  ##########################################################################################
  
  if (is.null(all_snps_wb)) {
    
    all_snps_wb <- data.frame(
      
      snp = character(0),
      
      chromosome = character(0),
      
      position = numeric(0),
      
      qvalue = numeric(0)
      
    )
    
  }
  
  
  ##########################################################################################
  # WRITE SNP TABLE
  ##########################################################################################
  
  print(
    paste(
      "Writing table for chemical:",
      chemical
    )
  )
  
  
  file_name <- paste0(
    "RANDOM_all_snps_wb_4km_",
    chemical,
    ".txt"
  )
  
  
  folder_path <- paste0(
    "~\\sensitivity_analysis_4km"
  )
  
  
  file_path <- file.path(
    folder_path,
    file_name
  )
  
  
  write.table(
    
    all_snps_wb,
    
    file = file_path,
    
    row.names = FALSE,
    
    quote = FALSE,
    
    sep = "\t"
    
  )
  
  
  print(
    paste(
      "Finished:",
      chemical
    )
  )
  
  
}   # END CHEMICAL LOOP


############################################################################################
######################################## END ###############################################
############################################################################################