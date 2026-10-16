library(openxlsx)

# Load data for Table S1

h1 <- read.delim("heterozygosity_individual.txt", header = TRUE, sep= " ") #FROM SCRIPT 6

# Load data for Table S2
h2 <- read.delim("heterozygosity_expected.txt", header = TRUE, sep= " ") #FROM SCRIPT 6
h3 <- read.delim("heterozygosity_observed.txt", header = TRUE, sep= " ") #FROM SCRIPT 6


# Load data for Table S3

vif_climate_included   <- read.csv("vif_results_climate_included.csv") #FROM SCRIPT 3

# Load data for Table S4

#vif_3                  <- read.csv("vif_results_climate_excluded_3km.csv") #FROM SCRIPT 3
vif_4                  <- read.csv("vif_results_climate_excluded_4km.csv") #FROM SCRIPT 3
#vif_5                  <- read.csv("vif_results_climate_excluded_5km.csv") #FROM SCRIPT 3

# Load data for Table S5
all_tox_wb_after_sens_duplicate_4km <-readRDS("all_tox_wb_after_sens_duplicate_4km.RData") #FROM SCRIPT 5

# Load data for Table S6
removed_links_4km <-readRDS("removed_links_4km.RData") #FROM SCRIPT 5

# Load data for Table S7
all_tox_wb_after_sens_duplicate_3km <-readRDS("all_tox_wb_after_sens_duplicate_3km.RData") #FROM SCRIPT 5

# Load data for Table S8
all_tox_wb_after_sens_duplicate_5km <-readRDS("all_tox_wb_after_sens_duplicate_5km.RData") #FROM SCRIPT 5

# Load data for Table S9
removed_links_3km <-readRDS("removed_links_3km.RData") #FROM SCRIPT 5

# Load data for Table S10
removed_links_5km <-readRDS("removed_links_5km.RData") #FROM SCRIPT 5

# Load data for Table S11
#constructed manually with script 7

# Load data for Table S12
si12   <- read.csv("revision_outputs/SITableS12.csv")

# Load data for Table S13
si13 <-  read.csv("revision_outputs/SITableS13_common_3km_4km_5km_enrichment.csv")

# Load data for Table S14
si14 <-  read.csv("revision_outputs/SITableS14.csv")

# Load data for Table S15
si15 <-  read.csv("revision_outputs/SITableS15_class_enrichment.csv")

# Load data for Table S16
si16 <-  read.csv("revision_outputs/SITableS16.csv")

# Load data for Table S17
si17 <-  read.csv("revision_outputs/TableS17_CTD_LG_overlap.csv")

#################################################################################################################
#################################################################################################################
#################################################################################################################

wb <- createWorkbook()

#Tables S1

addWorksheet(wb, "Table S1")
writeData(wb, "Table S1", h1, keepNA = TRUE, na.string = "NA" )

#Table S2

addWorksheet(wb, "Table S2 exp")
writeData(wb, "Table S2 exp", h2, keepNA = TRUE, na.string = "NA" )

addWorksheet(wb, "Table S2 obs")
writeData(wb, "Table S2 obs", h3, keepNA = TRUE, na.string = "NA" )

#Table S3

addWorksheet(wb, "Table S3")
writeData(wb, "Table S3", vif_climate_included)

#Table S4

addWorksheet(wb, "Table S4")
writeData(wb, "Table S4", vif_4)

#Table S5

addWorksheet(wb, "Table S5")
writeData(wb, "Table S5", all_tox_wb_after_sens_duplicate_4km, keepNA = TRUE, na.string = "NA" )

#Table S6

addWorksheet(wb, "Table S6")
writeData(wb, "Table S6", removed_links_4km, keepNA = TRUE, na.string = "NA" )

#Table S7

addWorksheet(wb, "Table S7")
writeData(wb, "Table S7", all_tox_wb_after_sens_duplicate_3km, keepNA = TRUE, na.string = "NA" )

#Table S8

addWorksheet(wb, "Table S8")
writeData(wb, "Table S8", all_tox_wb_after_sens_duplicate_5km, keepNA = TRUE, na.string = "NA" )

#Table S9

addWorksheet(wb, "Table S9")
writeData(wb, "Table S9", removed_links_3km, keepNA = TRUE, na.string = "NA" )

#Table S10

addWorksheet(wb, "Table S10")
writeData(wb, "Table S10", removed_links_5km, keepNA = TRUE, na.string = "NA" )

#Table S11


#Table S12

addWorksheet(wb, "Table S12")
writeData(wb, "Table S12", si12, keepNA = TRUE, na.string = "NA" )

#Table S13

addWorksheet(wb, "Table S13")
writeData(wb, "Table S13", si13, keepNA = TRUE, na.string = "NA" )

#Table S14

addWorksheet(wb, "Table S14")
writeData(wb, "Table S14", si14, keepNA = TRUE, na.string = "NA" )

#Table S15

addWorksheet(wb, "Table S15")
writeData(wb, "Table S15", si15, keepNA = TRUE, na.string = "NA" )

#Table S16

addWorksheet(wb, "Table S16")
writeData(wb, "Table S16", si16, keepNA = TRUE, na.string = "NA" )

#Table S17

addWorksheet(wb, "Table S17")
writeData(wb, "Table S17", si17, keepNA = TRUE, na.string = "NA" )

#################################################################################################################
#################################################################################################################
#################################################################################################################

#FINAL EXCEL

saveWorkbook(
  wb,
  "SI Tables S1 - S17.xlsx",
  overwrite = TRUE
)
