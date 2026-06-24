# ==============================================================================
# basemaps.R
# ==============================================================================
# Purpose: Loads shapefiles and defines basemap plotting functions for 
#          USUV phylodynamics in the Netherlands visualizations 
# Usage: Source this file from data visualization script to access basemaps.
#
# Shapefile sources:
#   - NL_ADM0_Land.shp, provincie_gegeneraliseerd.shp:
#       CBS, Kadaster, "CBSGebiedsindelingen2023"
#       https://service.pdok.nl/cbs/gebiedsindelingen/atom/v1_0/index.xml
#   - CNTR_RG_01M_2020_4326_NLAdapt.shp:
#       EuroGeographics - "Countries" (01M, 2020); NL boundary replaced with the CBS/Kadaster boundary above for exact alignment.
#       https://ec.europa.eu/eurostat/web/gisco/geodata/administrative-units/countries
#   - NLRiversEuMerged.shp 
#       (rivers, Netherlands portion):
#       Rijkswaterstaat (2024) "Kaderrichtlijn Water Oppervlaktewaterlichamen (vlakken)"
#       https://maps.rijkswaterstaat.nl/dataregister-publicatie/srv/fre/catalog.search#/metadata/rws1680f-68b5-4ff3-94a4-9c24109ffd5e
#       (rivers, neighbouring countries):
#       European Environment Agency (2020), EU-Hydro River Network Database
#       2006-2012 (vector), Europe - version 1.3, Copernicus Land Monitoring Service
#       https://doi.org/10.2909/393359a7-7ebd-4a52-80ac-1a18d5f3db9c
# ==============================================================================
# Set base_dir to your project root before running
# (if this file is sourced from another script that already set base_dir,
# this will just re-confirm the same working directory)
base_dir <- "."
setwd(base_dir)

# Load required libraries
library(terra)
library(sf)
library(foreign)
library(units)

# ==============================================================================
# 1. DEFINE SPATIAL EXTENTS
# ==============================================================================

## 1.1 Netherlands Extent (with buffer for neighboring countries)
# ------------------------------------------------------------------------------
# Load Netherlands administrative boundary (country level)
NLadm0 <- vect("Shapefiles/NL_ADM0_Land.shp")

# Get bounding box and expand it by 0.1 degrees in all directions
NLbbox <- st_bbox(NLadm0)
exp_NLbbox <- NLbbox

# Add buffer to include neighboring regions
exp_NLbbox["xmin"] <- NLbbox["xmin"] - 0.1 
exp_NLbbox["xmax"] <- NLbbox["xmax"] + 0.1 
exp_NLbbox["ymin"] <- NLbbox["ymin"] - 0.1 
exp_NLbbox["ymax"] <- NLbbox["ymax"] + 0.1

# Convert expanded bounding box to spatial object
NLextent <- st_as_sfc(st_bbox(exp_NLbbox), crs = 4326)

## 1.2 European Extent (for broader regional context)
# ------------------------------------------------------------------------------
# Define larger extent covering Western Europe
EUbbox <- st_bbox(
  c(xmin = -12.2, ymin = 35, xmax = 27, ymax = 61.2), 
  crs = st_crs(4326)
)
EUextent <- st_as_sfc(EUbbox, crs = 4326)

# ==============================================================================
# 2. LOAD AND PREPARE SHAPEFILES
# ==============================================================================

## 2.1 Netherlands Provincial Boundaries (ADM1)
# ------------------------------------------------------------------------------
# Load provincial boundaries of the Netherlands
NLadm1 <- vect("Shapefiles/provincie_gegeneraliseerd.shp")

# Add two-letter province codes from lookup table
ProvCode <- read.csv("Shapefiles/TwoLettersProvinceCodeNL.csv", header = TRUE, sep = ";", stringsAsFactors = FALSE)

# Extract unique province codes and merge with shapefile
ProvCode <- unique(ProvCode[, c("statcode", "X2Lcode")])
NLadm1 <- terra::merge(NLadm1, ProvCode, by.x = "statcode", by.y = "statcode", all.x = TRUE)

## 2.2 European Country Boundaries
# ------------------------------------------------------------------------------
# Load European administrative boundaries
EUadm0 <- vect("Shapefiles/CNTR_RG_01M_2020_4326_NLAdapt.shp")

# Crop to Netherlands extent (for detailed regional maps)
EUadm0NL <- crop(EUadm0, exp_NLbbox)

# Crop to broader European extent (for continental context)
EUadm0.1 <- crop(EUadm0, EUbbox)

## 2.3 Rivers
# ------------------------------------------------------------------------------
# Load river shapefiles and crop to Netherlands extent
Riv <- vect("Shapefiles/NLRiversEuMerged.shp")
Riv <- crop(Riv, exp_NLbbox)

# ==============================================================================
# 3. BASEMAP PLOTTING FUNCTIONS
# ==============================================================================

## 3.1 Netherlands with Rivers and Neighboring Countries
# ------------------------------------------------------------------------------
BaseMap.Riv <- expression({
  # Set plot parameters: no margins, no axis padding
  par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  
  # Plot base extent with water-themed background
  plot(NLextent, col = "#eef5f7", border = "#737373", lwd = 0.4)
  # Add neighboring countries 
  plot(EUadm0NL, col = "#F2F2F2", border = NA, lwd = 0.1, add = TRUE)
  # Add Netherlands land
  plot(NLadm0, col = "#e9e0de", border = NA, lwd = 0.1, add = TRUE)
  # Add rivers
  plot(Riv, col = "#cbd5d8", border = "#acbdc1", lwd = 0.1, add = TRUE)
  # Add country borders
  plot(EUadm0NL, col = NA, border = "#737373", lwd = 0.4, add = TRUE)
})

#eval(BaseMap.Riv) # to check

## 3.2 Netherlands Only
# ------------------------------------------------------------------------------
BaseMap <- expression({
  # Set plot parameters: no margins, no padding, no bounding box
  par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i", bty = "n")
  
  # Plot Netherlands
  plot(NLadm0, col = "#e9e0de", border = NA, axes = FALSE,box=FALSE)

  # Add provincial boundaries 
  plot(NLadm1, col = NA, border = "#bdbdbd", lwd = 0.3, add = TRUE)
  # Add national border 
  plot(NLadm0, col = NA, border = "#737373", lwd = 0.4, add = TRUE)
})

#eval(BaseMap) # to check

## 3.3 Netherlands Only with Provincial Labels
# ------------------------------------------------------------------------------
BaseMapLab <- expression({
  # Set plot parameters: no margins, no padding, no bounding box
  par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i", bty = "n")
  
  # Plot Netherlands
  plot(NLadm0, col = "#e9e0de", border = NA, lwd = 0.1, axes = FALSE)
  # Add provincial boundaries 
  plot(NLadm1, col = NA, border = "#bdbdbd", lwd = 0.3, add = TRUE)
  # Add national border 
  plot(NLadm0, col = NA, border = "#737373", lwd = 0.4, add = TRUE)
  # Add two-letter province codes as labels
  text(NLadm1, labels = NLadm1$X2Lcode, col = "#737373", cex = 1, pos = 4, offset = 0.05
  )
})

#eval(BaseMapLab) # to check

## 3.4 Northwestern Europe with Netherlands Highlighted
# ------------------------------------------------------------------------------
EUMap <- expression({
  # Set plot parameters: no margins, no padding
  par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  
  # Plot European extent with white background
  plot(EUextent, col = "#FFFFFF", border = "#737373")
  # Add European countries
  plot(EUadm0.1, col = "#F2F2F2", border = NA, lwd = 0.1, add = TRUE)
  # Highlight Netherlands 
  plot(NLadm0, col = "#e9e0de", border = NA, lwd = 0.1, add = TRUE)
  # Add country borders
  plot(EUadm0.1, col = NA, border = "#737373", lwd = 0.4, add = TRUE)
})

#eval(EUMap) # to check

# ==============================================================================
# END OF SCRIPT
# ==============================================================================