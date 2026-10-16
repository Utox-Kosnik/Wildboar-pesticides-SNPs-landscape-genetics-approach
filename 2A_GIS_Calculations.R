Sys.setenv(LANG = "en")

rm(list = ls())
library(sf)
library(terra)
library(exactextractr)
library(data.table)

#######################################################################################################################################
################################                        ################################################################################
################################        CATCHMENT       #####################################################################################################
################################                        #######################################################################################################
#######################################################################################################################################


# ------------------------------------------------------------
# 1. Load catchment raster
# ------------------------------------------------------------

catchments <- rast("catchments1.tif")
crs(catchments)
# Should be EPSG:3035

metadata_simple_wb <- read.table("metadata_simple_wb.txt", header = T) #FROM WB GENETIC DATA


# ------------------------------------------------------------
# 2. Convert coordinate dataframe to sf
# ------------------------------------------------------------

points <- st_as_sf(
  metadata_simple_wb,
  coords = c("longitude", "latitude"),
  crs = 4326,        # original longitude / latitude
  remove = FALSE
)

# Transform to ETRS89 / LAEA Europe
points_3035 <- st_transform(points, 3035)

# ------------------------------------------------------------
# 3. Create 4-km habitat circles
# ------------------------------------------------------------

habitat <- st_buffer(
  points_3035,
  dist = 4000          # EPSG:3035 is in metres
)

# ------------------------------------------------------------
# 4. Extract raster cells intersecting each habitat circle
# ------------------------------------------------------------

x <- exact_extract(
  catchments,
  habitat,
  include_cols = c("nr", "id", "pop_name", "country")
)

head(x[[1]])
names(x[[1]])

value_col <- setdiff(
  names(x[[1]]),
  c(
    "coverage_fraction",
    "nr",
    "id",
    "pop_name",
    "country"
  )
)

value_col

# ------------------------------------------------------------
# 5. Calculate area per catchment within each habitat circle
# ------------------------------------------------------------

cell_area_m2 <- prod(res(catchments)) #100 x 100 grid
cell_area_m2

catchment_area <- rbindlist(
  lapply(seq_along(x), function(i) {
    
    d <- as.data.table(x[[i]])
    
    # Rename raster-value column to catchment_ID
    setnames(d, "value", "catchment_ID")
    
    # Area of each raster cell actually inside the circle
    d[, area_m2 := coverage_fraction * cell_area_m2]
    
    # Sum area belonging to each catchment
    d[
      !is.na(catchment_ID),
      .(
        area_m2 = sum(area_m2),
        area_km2 = sum(area_m2) / 1e6
      ),
      by = .(nr,
             id,
             pop_name,
             country,
             catchment_ID)
    ]
  })
)
head(catchment_area)
#quality check, should be approx # 50.27 km2
area_check <- catchment_area[
  ,
  .(total_area_km2 = sum(area_km2)),
  by = .(nr, id)
]

summary(area_check$total_area_km2)
fwrite(catchment_area, "catchment_area_4km.csv")

# Current R buffers
as.numeric(st_area(habitat[1, ])) / 1e6

#######################################################################################################################################
################################                        ################################################################################
################################       LANDCOVER        #####################################################################################################
################################                        #######################################################################################################
#######################################################################################################################################

rm(list = ls())

# ------------------------------------------------------------
# 1. Load land use raster
# ------------------------------------------------------------

landuse <- rast("2018_landcover.tif")
crs(landuse)
# Should be EPSG:3035

metadata_simple_wb <- read.table("metadata_simple_wb.txt", header = T) #FROM WB GENETIC DATA

# ------------------------------------------------------------
# 2. Convert coordinate dataframe to sf
# ------------------------------------------------------------

points <- st_as_sf(
  metadata_simple_wb,
  coords = c("longitude", "latitude"),
  crs = 4326,        # original longitude / latitude
  remove = FALSE
)

# Transform to ETRS89 / LAEA Europe
points_3035 <- st_transform(points, 3035)

# ------------------------------------------------------------
# 3. Create 4-km habitat circles
# ------------------------------------------------------------

habitat <- st_buffer(
  points_3035,
  dist = 4000          # EPSG:3035 is in metres
)

# ------------------------------------------------------------
# 4. Extract raster cells intersecting each habitat circle
# ------------------------------------------------------------

y <- exact_extract(
  landuse,
  habitat,
  include_cols = c("nr", "id", "pop_name", "country")
)

head(y[[1]])
names(y[[1]])

value_col_y <- setdiff(
  names(y[[1]]),
  c(
    "coverage_fraction",
    "nr",
    "id",
    "pop_name",
    "country"
  )
)

value_col_y

# ------------------------------------------------------------
# 5. Calculate area per catchment within each habitat circle
# ------------------------------------------------------------

landuse_cell_area_m2 <- prod(res(landuse)) #100 x 100 grid
landuse_cell_area_m2

landuse_area <- rbindlist(
  lapply(seq_along(y), function(i) {
    
    d <- as.data.table(y[[i]])
    
    # Rename raster-value column to catchment_ID
    setnames(d, "value", "catchment_ID")
    
    # Area of each raster cell actually inside the circle
    d[, area_m2 := coverage_fraction * landuse_cell_area_m2]
    
    # Sum area belonging to each catchment
    d[
      !is.na(catchment_ID),
      .(
        area_m2 = sum(area_m2),
        area_km2 = sum(area_m2) / 1e6
      ),
      by = .(nr,
             id,
             pop_name,
             country,
             catchment_ID)
    ]
  })
)
head(landuse_area)
#quality check, should be approx # 50.27 km2 (4km radius), 78.5 km2 (5km radius)
landuse_area_check <- landuse_area[
  ,
  .(total_area_km2 = sum(area_km2)),
  by = .(nr, id)
]

summary(landuse_area_check$total_area_km2)


landuse_percentage <- landuse_area[, .(
  total_area_km2 = sum(area_km2, na.rm = TRUE),
  agricultural_area_km2 = sum(
    area_km2[catchment_ID %in% 12:22],
    na.rm = TRUE
  )
), by = .(nr, id, pop_name, country)]

landuse_percentage[, agricultural_percentage :=
                     agricultural_area_km2 / total_area_km2 * 100]

fwrite(landuse_percentage, "landuse_percentage_4km.csv")

#######################################################################################################################################
################################                        ################################################################################
################################    HUMAN POPULATION    #####################################################################################################
################################                        #######################################################################################################
#######################################################################################################################################

rm(list = ls())

library(sf)
library(terra)
library(exactextractr)
library(data.table)

# ------------------------------------------------------------
# 1. Load land use raster
# ------------------------------------------------------------

population <- rast("human_population_data\\JRC_1K_POP_2018.tif")
crs(population)
res(population)
# Should be EPSG:3035

metadata_simple_wb <- read.table("metadata_simple_wb.txt", header = T) #FROM WB GENETIC DATA

# ------------------------------------------------------------
# 2. Convert coordinate dataframe to sf
# ------------------------------------------------------------

points <- st_as_sf(
  metadata_simple_wb,
  coords = c("longitude", "latitude"),
  crs = 4326,        # original longitude / latitude
  remove = FALSE
)

# Transform to ETRS89 / LAEA Europe
points_3035 <- st_transform(points, 3035)

# ------------------------------------------------------------
# 3. Function to calculate population within habitat circle
# ------------------------------------------------------------

calculate_population <- function(radius_m) {
  
  # Create habitat circles
  habitat <- st_buffer(
    points_3035,
    dist = radius_m
  )
  
  # Extract population raster cells intersecting habitat circles
  pop_extract <- exact_extract(
    population,
    habitat,
    include_cols = c("nr", "id", "pop_name", "country")
  )
  
  # Combine results for all individuals
  population_total <- rbindlist(
    lapply(seq_along(pop_extract), function(i) {
      
      d <- as.data.table(pop_extract[[i]])
      
      # Rename raster value column
      setnames(d, "value", "population")
      
      # Population contribution of cell within habitat circle
      d[, population_inside := population * coverage_fraction]
      
      # Sum population within habitat circle
      d[
        !is.na(population),
        .(
          human_population = sum(population_inside, na.rm = TRUE)
        ),
        by = .(
          nr,
          id,
          pop_name,
          country
        )
      ]
    })
  )
  
  return(population_total)
}

# ------------------------------------------------------------
# 4. Calculate population for 3-, 4- and 5-km radius
# ------------------------------------------------------------

population_3km <- calculate_population(3000)
population_4km <- calculate_population(4000)
population_5km <- calculate_population(5000)


# ------------------------------------------------------------
# 5. Inspect results
# ------------------------------------------------------------

head(population_3km)

summary(population_3km$human_population)
summary(population_4km$human_population)
summary(population_5km$human_population)


# ------------------------------------------------------------
# 6. Save
# ------------------------------------------------------------

fwrite(population_3km, "population_3km.csv")
fwrite(population_4km, "population_4km.csv")
fwrite(population_5km, "population_5km.csv")

#######################################################################################################################################
################################                        ################################################################################
################################    CLIMATE ZONES       #####################################################################################################
################################                        #######################################################################################################
#######################################################################################################################################

rm(list = ls())

library(sf)
library(terra)
library(exactextractr)
library(data.table)

# ------------------------------------------------------------
# 1. Load climate raster
# ------------------------------------------------------------

climate <- rast("climate_data_1991_2020\\koppen_geiger_1p0_111km.tif")
crs(climate)
res(climate)
names(climate) <- "climate_class"

# ------------------------------------------------------------
# 2. Climate class lookup table
# ------------------------------------------------------------

climate_lookup <- data.table(
  climate_class = 1:30,
  climate_zone = c("Af", "Am", "Aw","BWh", "BWk", "BSh", "BSk","Csa", "Csb", "Csc","Cwa", "Cwb", "Cwc",
    "Cfa", "Cfb", "Cfc","Dsa", "Dsb", "Dsc", "Dsd","Dwa", "Dwb", "Dwc", "Dwd","Dfa", "Dfb", "Dfc", "Dfd","ET", "EF"))

# ------------------------------------------------------------
# 3. Load sampling coordinates
# ------------------------------------------------------------

metadata_simple_wb <- read.table("metadata_simple_wb.txt",header = TRUE)
points <- st_as_sf(metadata_simple_wb,coords = c("longitude", "latitude"),crs = 4326,remove = FALSE)

# Make buffers in metres
points_3035 <- st_transform(points, 3035)

# ------------------------------------------------------------
# 4. Raster-cell areas
# ------------------------------------------------------------

# Climate raster is WGS84, so cells are defined in degrees.
# Calculate their actual areas in km2.

climate_cell_area <- cellSize(climate,unit = "km")
names(climate_cell_area) <- "cell_area_km2"
climate_stack <- c(climate,climate_cell_area)

# ------------------------------------------------------------
# 5. Function to calculate climate zones within habitat circle
# ------------------------------------------------------------

calculate_climate <- function(radius_m) {
  
  # Create habitat circles
  habitat <- st_buffer(
    points_3035,
    dist = radius_m
  )
  
  # Transform circles to climate raster CRS
  habitat_climate <- st_transform(
    habitat,
    crs(climate)
  )
  
  # Extract intersecting raster cells
  z <- exact_extract(
    climate_stack,
    habitat_climate,
    include_cols = c(
      "nr",
      "id",
      "pop_name",
      "country"
    )
  )
  
  
  # ----------------------------------------------------------
  # Area of each climate class within circle
  # ----------------------------------------------------------
  
  climate_area <- rbindlist(
    
    lapply(seq_along(z), function(i) {
      
      d <- as.data.table(z[[i]])
      
      # Area of raster cell actually inside habitat circle
      d[, area_km2 :=
          coverage_fraction * cell_area_km2]
      
      # Sum area for each climate class
      d[
        !is.na(climate_class),
        .(
          area_km2 = sum(area_km2, na.rm = TRUE)
        ),
        by = .(nr,id,pop_name,country,climate_class)]
    })
  )
  
  
  # ----------------------------------------------------------
  # Add climate names
  # ----------------------------------------------------------
  
  climate_area <- climate_lookup[climate_area,on = "climate_class"]
  
  # ----------------------------------------------------------
  # Percentage of habitat circle per climate class
  # ----------------------------------------------------------
  
  climate_area[,
    climate_percentage :=area_km2 / sum(area_km2) * 100,
    by = .(nr,id,pop_name,country)]
  
  
  # ----------------------------------------------------------
  # Summary per wild boar
  # ----------------------------------------------------------
  
  climate_summary <- climate_area[
    ,
    .(
      # How many climate zones occur?
      n_climate_types = uniqueN(climate_zone),
      
      # Which climate zones occur?
      climate_types = paste(
        climate_zone[order(-area_km2)],
        collapse = ", "
      ),
      
      # Main climate zone
      dominant_climate =
        climate_zone[which.max(area_km2)],
      
      # Percentage occupied by main zone
      dominant_climate_percentage =
        climate_percentage[which.max(area_km2)]
      
    ),
    by = .(
      nr,
      id,
      pop_name,
      country
    )
  ]
  
  
  return(
    list(
      summary = climate_summary,
      area = climate_area
    )
  )
}

# ------------------------------------------------------------
# 4. Calculate population for 10-km radius
# ------------------------------------------------------------

climate_10km <- calculate_climate(10000)

# ------------------------------------------------------------
# 5. Inspect results
# ------------------------------------------------------------

table(climate_10km$summary$n_climate_types) #35 individuals have 2 climate zones

multiple_10km <- climate_10km$summary[n_climate_types > 1,.(nr, id)]
climate_10km$area[multiple_10km,on = .(nr, id)][order(nr, -climate_percentage)]

#34 out of 35 individuals have a clear major (>75%) vs minor climate zone, except 1 individual that has 56 vs 44%

climate_10km_major <- climate_10km$summary[
  , .(nr,id,pop_name,country,climate_zone = dominant_climate,climate_percentage = dominant_climate_percentage)]

# ------------------------------------------------------------
# 6. Save
# ------------------------------------------------------------

fwrite(climate_10km_major, "climate_10km.csv")

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
rm(list = ls())

catchment_area_3km <- fread("catchment_area_3km.csv")
catchment_area_4km <- fread("catchment_area_4km.csv")
catchment_area_5km <- fread("catchment_area_5km.csv")

landuse_percentage_3km <- fread("landuse_percentage_3km.csv")
landuse_percentage_4km <- fread("landuse_percentage_4km.csv")
landuse_percentage_5km <- fread("landuse_percentage_5km.csv")

population_3km <- fread("population_3km.csv")
population_4km <- fread("population_4km.csv")
population_5km <- fread("population_5km.csv")

climate_10km <- fread("climate_10km.csv")

# ------------------------------------------------------------
# 3 km
# ------------------------------------------------------------

data_3km <- catchment_area_3km[
  landuse_percentage_3km,
  agricultural_percentage := i.agricultural_percentage,
  on = .(nr, id, pop_name, country)
]

data_3km <- data_3km[
  population_3km,
  human_population := i.human_population,
  on = .(nr, id, pop_name, country)
]

data_3km <- data_3km[
  climate_10km,
  `:=`(
    climate_zone = i.climate_zone,
    climate_percentage = i.climate_percentage
  ),
  on = .(nr, id, pop_name, country)
]



# ------------------------------------------------------------
# 4 km
# ------------------------------------------------------------

data_4km <- catchment_area_4km[
  landuse_percentage_4km,
  agricultural_percentage := i.agricultural_percentage,
  on = .(nr, id, pop_name, country)
]

data_4km <- data_4km[
  population_4km,
  human_population := i.human_population,
  on = .(nr, id, pop_name, country)
]

data_4km <- data_4km[
  climate_10km,
  `:=`(
    climate_zone = i.climate_zone,
    climate_percentage = i.climate_percentage
  ),
  on = .(nr, id, pop_name, country)
]


# ------------------------------------------------------------
# 5 km
# ------------------------------------------------------------

data_5km <- catchment_area_5km[
  landuse_percentage_5km,
  agricultural_percentage := i.agricultural_percentage,
  on = .(nr, id, pop_name, country)
]

data_5km <- data_5km[
  population_5km,
  human_population := i.human_population,
  on = .(nr, id, pop_name, country)
]

data_5km <- data_5km[
  climate_10km,
  `:=`(
    climate_zone = i.climate_zone,
    climate_percentage = i.climate_percentage
  ),
  on = .(nr, id, pop_name, country)
]


# Save
fwrite(data_3km, "data_3km.csv")
fwrite(data_4km, "data_4km.csv")
fwrite(data_5km, "data_5km.csv")
