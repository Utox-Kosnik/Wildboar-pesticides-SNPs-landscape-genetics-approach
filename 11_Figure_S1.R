library(data.table)
library(ggplot2)
library(patchwork)

packageVersion("patchwork")

rm(list = ls())

#PLOT 1 : Sum up all the values across all countries for each year and pesticide category

# ------------------------------------------------------------
# 1. Read FAOSTAT data
# ------------------------------------------------------------

file <- "FAOSTAT_data_en_9-16-2026.csv"

# FAOSTAT file has unusual outer quotation marks
raw <- readLines(file)

header <- raw[1]

dat_lines <- raw[-1]

# Remove outer quotes and replace doubled quotes
dat_lines <- sub('^"', '', dat_lines)
dat_lines <- sub('"$', '', dat_lines)
dat_lines <- gsub('""', '"', dat_lines, fixed = TRUE)

# Read cleaned data
dat <- fread(
  text = paste(c(header, dat_lines), collapse = "\n")
)

# ------------------------------------------------------------
# 2. Convert variables
# ------------------------------------------------------------

dat[, Year := as.integer(Year)]
dat[, Value := as.numeric(Value)]

# Check the four categories
unique(dat$Item)

# ------------------------------------------------------------
# 3. Sum pesticide use across countries per year
# ------------------------------------------------------------

plot_dat <- dat[
  ,
  .(
    pesticide_use = sum(Value, na.rm = TRUE)
  ),
  by = .(Year, Item)
]

# ------------------------------------------------------------
# 4. Plot
# ------------------------------------------------------------

p <- ggplot(
  plot_dat,
  aes(
    x = Year,
    y = pesticide_use,
    colour = Item
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  labs(
    x = "Year",
    y = "Pesticide use (tonnes)",
    colour = "Pesticide category"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  )

p

#PLOT 2 : Country specific trends

# ------------------------------------------------------------
# Prepare data
# ------------------------------------------------------------

dat[, Year := as.integer(Year)]
dat[, Value := as.numeric(Value)]

# Check countries and categories
unique(dat$Area)
unique(dat$Item)

# ------------------------------------------------------------
# Plot: individual country trends
# ------------------------------------------------------------

p_country <- ggplot(
  dat,
  aes(
    x = Year,
    y = Value,
    group = Area,
    colour = Area
  )
) +
  geom_line(
    linewidth = 0.7,
    alpha = 0.8
  ) +
  
  facet_wrap(
    ~ Item,
    ncol = 2,
    scales = "free_y"
  ) +
  
  scale_y_log10() +
  
  labs(
    x = "Year",
    y = "Pesticide use (tonnes)",
    colour = "Country"
  ) +
  
  theme_bw() +
  
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 11),
    
    strip.text = element_text(
      size = 13,
      face = "bold"
    ),
    
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9),
    
    panel.grid.minor = element_blank()
  )

p_country



combined_plot <- p / p_country +
  plot_layout(
    heights = c(1, 2)
  )

combined_plot

