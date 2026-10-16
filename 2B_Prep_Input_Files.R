############################################################################
#Prepping input files for RDA: chemical list + matrix + metadata
#based on 4 km home ranges 
###########################################################################

rm(list = ls()) 

######################
#Packages needed 
######################

library(data.table)
library(dplyr) #%>%
library(tidyr) #pivot_wider

######################
#Input files needed
######################

wb_matrix_no_na <- readRDS("final_matrix_wb_europe_no_na.RData") #FROM SCRIPT 1
#306 individuals, 40700 SNPs
colnames(wb_matrix_no_na) <- c(1:40700)

tox_names<-read.table("tox.txt",sep = "\t", header = T) #FROM PISTOCCHI

#qgis <-read.table("qgis_pisto_wb_4km_radius_final.txt", header = T) #FROM QGIS
data_4km <- fread("data_4km.csv") #FROM PREVIOUS R SCRIPT (either data_4km, data_5km or data_3km)
#CATCHMENT_ID = polygon identity number (each polygon has a known concentration) + area_m2 = how much m2 of this polygon was inside the 4km radius

pest_conc<-read.table("calc_conc.txt",sep = "\t", header = T) #FROM PISTOCCHI

metadata_simple_wb <- read.table("metadata_simple_wb.txt", header = T) #FROM WB GENETIC DATA

test_match_order2 <- function(x,y) {
  if (isTRUE(all.equal(x,y))) print('Perfect match in same order')
  if (!isTRUE(all.equal(x,y)) && isTRUE(all.equal(sort(x),sort(y)))) print('Perfect match in wrong order')
  if (!isTRUE(all.equal(x,y)) && !isTRUE(all.equal(sort(x),sort(y)))) print('No match')
}

###########################################################################

######################
# PREP CHEMICAL NAMES (long name + Xcode) --> chemical list
######################

colnames(tox_names) <- c("xtox_name","name","extra")
tox_names <- tox_names %>% mutate(name = recode(name,
                                                "1-Naphthylacetic acid (1-NAA)" = "1-Naphthylacetic_acid",
                                                "Boscalid (formerly nicobifen)" = "Boscalid",
                                                "Chloridazon (aka pyrazone)" = "Chloridazon",
                                                "Diquat  (dibromide)" = "Diquat",
                                                "Flufenacet (formerly fluthiamide)" = "Flufenacet",
                                                "Glyphosate (incl trimesium aka sulfosate)" = "Glyphosate",
                                                "Imazalil (aka enilconazole)" = "Imazalil",
                                                "Methiocarb (aka mercaptodimethur)" = "Methiocarb",
                                                "Rimsulfuron (aka renriduron)" = "Rimsulfuron"
))
tox_names[, 1] <- paste0("X", tox_names[, 1])
chemical_list <- tox_names$name #148 chemicals

#look up the 148 pesticides and the year of their first application
pest_inclusion <- read.csv("260426-FullPesticides_Annotation.csv") 
pest_yes <- pest_inclusion %>%
  filter(Approximated_year <= 2000)
pest_name_included_before_2000 <- pest_yes$tox_name #121 out of 148 were in use before 2000
chemical_list <- pest_name_included_before_2000 #121
chemical_list[1] <- "1-Naphthylacetic_acid" #just add _ character
saveRDS(chemical_list, "chemical_list.RData") #121

######################
#PREP CATCHMENT DATA (AREA AND RASTERVALUE PER INDIVIDUAL)
######################

data_4km_wide <- data_4km %>%
  select(id,catchment_ID,area_m2) %>%
  mutate(catchment_ID =as.character(catchment_ID),area_m2 =as.numeric(area_m2)) %>%
  pivot_wider(id_cols = id,names_from =catchment_ID, values_from =area_m2,values_fn =sum, values_fill = 0)

# IDs
id_catchment <- data_4km_wide$id


# Remove ID column:
# rows = individuals
# columns = catchment IDs
catchment_area_values <-data_4km_wide[, -1]
catchment_area_values <-as.data.frame(catchment_area_values)

# Total mapped catchment area
# for each individual
row_sums <-rowSums(catchment_area_values,na.rm = TRUE)

# Check area in km2
summary(row_sums / 1e6)

# A complete 4 km radius circle should be approximately:
# pi * 4^2 = 50.27 km2


# Catchment IDs
# that need to be matched
# with pesticide concentrations
order_rastervalues <-colnames(catchment_area_values)

############################################################################
# HUMAN POPULATION
############################################################################

# Keep only one value per individual.
human_population <- unique(data_4km[,.(id,human_population)])

# Check that there is exactly one row per individual
nrow(human_population)
# Should be 282

# Check that IDs are unique
if (anyDuplicated(human_population$id)) {stop("Some individuals have more than one human_population value.")}

############################################################################
# AGRICULTURAL LAND PERCENTAGE
############################################################################

# agricultural_percentage is repeated
# for every catchment row of an individual.
#
# Keep only one value per individual.

agriculture <- unique(data_4km[,.(id,agricultural_percentage)])

# Check that there is exactly one row per individual
nrow(agriculture)

# Should be 282


# Check that IDs are unique
if (anyDuplicated(agriculture$id)) {stop("Some individuals have more than one agricultural_percentage value.")}

############################################################################
# CLIMATE
############################################################################

climate <- unique(data_4km[,.(id,climate_zone)])

# Check that there is exactly one row per individual
nrow(climate)
# Should be 282

######################
#PREP CONCENTRATION DATA (LINK BETWEEN RASTERVALUE AND CONCENTRATION AND Xcode)
######################

pest_conc$catch_name <- as.character(pest_conc$WSO1_ID)
help_file <- pest_conc[,c(names(pest_conc)[9:156],"catch_name")] #catch_name = raster value, 149 columns = 148 chemicals + 1 polygon identity column
help_file_filtered <- help_file[help_file$catch_name %in% order_rastervalues, ]  # Keep only the 1129 polygons that we need
help_file_ordered <- help_file_filtered[match(order_rastervalues, help_file_filtered$catch_name), ]
#columns are chemicals, rows are polygon identities

######################
# CHECK CATCHMENT MATCHING
######################

# There should not be any missing matches
sum(is.na(help_file_ordered$catch_name))

# Ideally returns:
# 0

# Check that order matches exactly
all(help_file_ordered$catch_name ==order_rastervalues)

# Should return:
# TRUE

######################
#MERGE CHEMICAL NAMME WITH QGIS DATA AND CONCENTRATION DATA
######################

#all chemical full names

all_conc_all_tox = NULL
for (chemical in chemical_list){
  print(chemical)
  #get the chemical Xcode
  help <- tox_names %>% filter(tox_names$name == chemical)
  xtox <- help$xtox_name 
  #FILTER THE COLUMNS WITH CONC AND RASTERVALUES FOR 1 CHEMICAL (Xcode)
  conc_values <- as.numeric(help_file_ordered[[xtox]])
  
  conc_final = NULL
  for (i in 1:(nrow(catchment_area_values))){ #iterate through rows 
    all_conc <- 0
    all_area <- 0
    for (j in 1:ncol(catchment_area_values)){ #iterate through columns 
      area_num <- as.numeric(catchment_area_values[i,j])
      conc_num <- conc_values[j]
      if ((area_num != 0)&& !is.na(conc_num)){
        conc <- area_num*conc_num
        all_conc <- all_conc + conc
        all_area <- all_area + area_num
      }else{
      }
    }
    id_conc = all_conc / all_area
    conc_final <- rbind(conc_final, id_conc)  
  }
  all_conc_all_tox <- cbind(all_conc_all_tox, conc_final)
  
}
colnames(all_conc_all_tox) <- chemical_list
metadata <-data.frame(id =id_catchment, all_conc_all_tox,check.names = FALSE)

final_metadata <- metadata_simple_wb %>%
  inner_join(metadata, by = "id") %>%
  left_join(as.data.frame(agriculture), by = "id") %>%
  left_join(as.data.frame(climate), by = "id")%>%
  left_join(as.data.frame(human_population), by = "id")

######################
# CHECK
######################

dim(final_metadata)

names(final_metadata)

summary(final_metadata$agricultural_percentage)
table(final_metadata$climate_zone)
######################
#ADJUST MATRIX
######################

ids_to_keep <- final_metadata$id

wb_matrix_no_na_282 <- wb_matrix_no_na[rownames(wb_matrix_no_na) %in% ids_to_keep, ] 
test_match_order2(rownames(wb_matrix_no_na_282), ids_to_keep) # "Perfect match in wrong order"
right_order <- rownames(wb_matrix_no_na_282)
wb_final_metadata <- final_metadata[match(right_order, final_metadata$id),]
test_match_order2(rownames(wb_matrix_no_na_282), (wb_final_metadata$id))

#replaces Inf, -Inf, and NaN values with NA
wb_final_metadata[sapply(wb_final_metadata, is.numeric)] <- lapply(wb_final_metadata[sapply(wb_final_metadata, is.numeric)], function(x) {
  x[is.infinite(x)] <- NA
  x[is.nan(x)] <- NA
  return(x)
})

str(wb_final_metadata)
wb_final_metadata$country<-as.factor(wb_final_metadata$country)
wb_final_metadata$pop_name<-as.factor(wb_final_metadata$pop_name)

############################################################################
# FINAL CHECKS
############################################################################

str(wb_final_metadata)
dim(wb_matrix_no_na_282)
dim(wb_final_metadata)

# Matrix and metadata
# must contain same number of individuals
nrow(wb_matrix_no_na_282) == nrow(wb_final_metadata)

# IDs must be in exactly same order
all(rownames(wb_matrix_no_na_282) ==wb_final_metadata$id)

######################
#SAVE DATA
######################

saveRDS(wb_final_metadata, "wb_final_metadata_4km.Rdata") 
saveRDS(wb_matrix_no_na_282, "wb_matrix_no_na_282.RData")
