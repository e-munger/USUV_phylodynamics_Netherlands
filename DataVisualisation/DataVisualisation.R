# =============================================================================
# Title: Data visualisation - Phylodynamics of USUV in the Netherlands
# Description: Script used to generate Figures 1-5 of the manuscript "Phylodynamic analyses reveal persistence 
# and extensive intermixing of Usutu virus in the Netherlands"
# Author: Emmanuelle Munger
# Several sections (notably continuous and country wide discrete phylogeography visualisation, Figs. 3-4) 
# adapt code originally written by Fabiana Gambaro and Simon Dellicour.
# Last edit: 2026-06-01
#
# Companion scripts (must be in base_dir): CustomThemes.R, basemaps.R
# Paths are relative to `base_dir`; update folder names to match your setup.
#
# Input data (BEAST MCC/posterior trees, GenBank metadata, shapefiles, environmental rasters) not included — update paths to your own files.
# =============================================================================

# Set base_dir to your project root before running
base_dir <- "."
setwd(base_dir)

# Path to Figures repository
FigRepo <- file.path(base_dir, "Figures")

# Custom ggplot theme and basemap definitions used throughout figures
source(file.path(base_dir, "CustomThemes.R"))
source(file.path(base_dir, "basemaps.R"))

## Fig. 1: Global DTA & sequence data availability ------------------------

library(tidyverse)
library(treeio)
library(ggtree)
library(ggnewscale)
library(seraphim)
library(diagram)
library(sf)
library(terra)
library(ape)
library(patchwork)
library(scales)

# (i) Vizualize MCC trees of datasets NFG + E and NFG + NS5

# MCC trees were obtained by running TreeAnnotator on the posterior trees
# from the discrete phylogeographic analysis (BEAST v1.10), with 10% burn-in.
# Read MCC trees
DirectoryGlobDTA <- file.path(base_dir, "1_TreeTimeDTA_E_NS5")

MCCglE <- read.beast(file.path(DirectoryGlobDTA, "Global_NFG_E/USUVGlobal_E430_treetimeDTA.MCC.tree"))
MCCglNS5 <- read.beast(file.path(DirectoryGlobDTA, "Global_NFG_NS5/USUVGlobal_NS5last_treetimeDTA_MCC.tree"))

# Read metadata curated from GenBank
MetaU <- read.csv("GlobalMetadata/250606_MinUSUVAllMetadata_Country.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE) %>%
  dplyr::select(- X)

# Prepare the labbeling of tree and set the color scale
# Check how many locations there are
unique(MetaU$collection_countryname) #25 countries 
table(MetaU$collection_countryname) #Check number of sequences per country
sort(table(MetaU$collection_countryname), decreasing = TRUE)

# Newly generated sequences from Dutch surveillance (2016-2023) are distinguished
# from previously published Dutch sequences
MetaNL <- read.csv("GenBank_ENA_Submissions/260225_ENA_Metadata_USUV_SequencesToBeSubmitted.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE) %>%
  dplyr::select(host_subject_id, Included_In)

MetaU <- left_join(MetaU, MetaNL, by = c("SeqID" = "host_subject_id"))

MetaU <- MetaU %>%
  mutate(collection_countryname =
           case_when(Included_In == "MungerE_Phylodynamics_2016_2023" ~ "NetherlandsNew",
                     TRUE ~ collection_countryname))

# Countries with <10 sequences are collapsed into broad geographic regions
# to simplify visualization 
country_mapping <- c(
  "NetherlandsNew" = "NetherlandsNew",
  "Netherlands" = "Netherlands",
  "Italy" = "Italy",
  "Germany" = "Germany",
  "Austria" = "Europe other",
  "Czechia" = "Europe other",
  "Poland" = "Europe other",
  "Hungary" = "Europe other",
  "Slovakia" = "Europe other",
  "Belgium" = "Europe other",
  "Senegal" = "Africa",
  "France" = "Europe other",
  "Israel" = "Middle East",
  "Spain" = "Europe other",
  "Croatia" = "Europe other",
  "Greece" = "Europe other",
  "Uganda" = "Africa",
  "UnitedKingdom" = "Europe other",
  "CentralAfricanRepublic" = "Africa",
  "Luxembourg" = "Europe other",
  "Portugal" = "Europe other",
  "Romania" = "Europe other",
  "Serbia" = "Europe other",
  "SouthAfrica" = "Africa",
  "Sweden" = "Europe other",
  "Switzerland" = "Europe other")
country_mapping

unique_mapped_countries <- unique(country_mapping)
unique_mapped_countries #7 locations

#Define the order of the countries and specific colors
country_order <- c("Germany",
                   "Italy",
                   "Europe other",
                   "Africa",       
                   "Middle East",
                   "Netherlands",
                   "NetherlandsNew")  

colors <- c(
  "Germany" =  "#5b82b4", 
  "Italy"="#b5decc",
  "Europe other"="#C19B8E",
  "Africa"= "#3b3b77",
  "Middle East"="#ffbd62",
  "Netherlands" = "#cc1337",
  "NetherlandsNew" = NULL) # NetherlandsNew excluded — drawn separately

tipsize <- c(
  "Germany" =  0.9, 
  "Italy"= 0.9,
  "Europe other"= 0.9,
  "Africa"= 0.9,
  "Middle East"= 0.9,
  "Netherlands" = 0.9,
  "NetherlandsNew" = NULL)

# Add country class to metadata
MetaU <- MetaU %>%
  mutate(CountryClass = recode(collection_countryname, !!!country_mapping))

# Make sure CountryClass is a factor with the specified levels (order of countries)
MetaU$CountryClass <- factor(MetaU$CountryClass, levels = country_order)

# MCC tree for NFG+E dataset; branches coloured by inferred ancestral location
# (NL vs other) from the preliminary two-state discrete phylogeographic analysis;
# mrsd set to the most recent sampling date in the combined global dataset

t1 <- ggtree(MCCglE, mrsd="2024-10-12", ladderize=T, right = T, aes(color = location), size = 0.32) +
  theme_tree2()+
  geom_rootedge(0.8, color = "#a7a3e2") + #plots a root
  scale_color_manual(values = c("NL" = "#cc1337", "other" = "#a7a3e2")) 

# Annotate tips with sampling-location colour class; newly generated Dutch sequences
# drawn as circles with outlines to visually distinguish them;
# Lineages Europe 3 and Africa 3 are labelled
t1.1 <- t1  %<+% MetaU + new_scale_color() +
  geom_tippoint(aes(color = CountryClass, size = CountryClass), alpha = 0.6) +
  scale_color_manual("CountryClass", values = colors) +
  scale_size_manual("CountryClass", values = tipsize) +
  geom_tippoint(data = . %>% 
                  filter (CountryClass == "NetherlandsNew" & isTip == TRUE),
                color = "#731331", # stroke color 
                fill = "#cc133799", # "#cc1337" with 60% opacity
                size = 1.5, 
                shape = 21, stroke = 0.3) +
  Arbo_theme_grid() +
  theme(legend.position = "none",
        axis.text.y = element_blank(),      
        axis.ticks.y = element_blank(),     
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank()) +
  ggtitle(NULL, subtitle = "Near full genomes + E") +   # Subtitle only 
  scale_x_continuous(breaks = seq(1920, 2025, by = 10),limits = c(1920, 2025)) + 
  geom_text(aes(label=ifelse(node==1283, "Europe 3", "")), 
            hjust=1.1, vjust=-0.5, size=2, color="#737373") +
  geom_text(aes(label=ifelse(node==908, "Africa 3", "")), 
            hjust=1.1, vjust=-0.5, size=2, color="#737373")

# MCC tree for NFG+NS5 dataset; plotted with identical settings
t2 <- ggtree(MCCglNS5, mrsd="2024-10-12", ladderize=T, right = T, aes(color = location), size = 0.32) +
  theme_tree2()+
  geom_rootedge(0.8, color = "#a7a3e2") + #plots a root
  scale_color_manual(values = c("NL" = "#cc1337", "other" = "#a7a3e2")) 

t2.2 <- t2  %<+% MetaU + new_scale_color() +
  geom_tippoint(aes(color = CountryClass, size = CountryClass), alpha = 0.6) +
  scale_color_manual("CountryClass", values = colors) +
  scale_size_manual("CountryClass", values = tipsize) +
  geom_tippoint(data = . %>% 
                  filter (CountryClass == "NetherlandsNew" & isTip == TRUE),
                color = "#731331", # stroke color 
                fill = "#cc133799", # "#cc1337" with 60% opacity
                size = 1.5, 
                shape = 21, stroke = 0.3) +
  Arbo_theme_grid() +
  theme(legend.position = "none",
        axis.text.y = element_blank(),      
        axis.ticks.y = element_blank(),     
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank()) +
  ggtitle(NULL, subtitle = "Near full genomes + NS5") +
  scale_x_continuous(breaks = seq(1920, 2025, by = 10), limits = c(1920, 2025)) +
  geom_text(aes(label=ifelse(node==1278, "Europe 3", "")), 
            hjust=1.1, vjust=-0.5, size=2, color="#737373") +
  geom_text(aes(label=ifelse(node==950, "Africa 3", "")), 
            hjust=1.3, vjust=0.6, size=2, color="#737373")

# Export landscape A5
pdf(paste0(FigRepo,"GlobalDTA/","Trees","rev.pdf"), 
    width = 8.27, height = 5.83)

t1.1 | t2.2

dev.off()

# (i) Visualize the number of sequences available per country in NFG + E and NFG + NS5 datasets 

# Prepare metadata for maps
# African and Middle Eastern sequences are plotted below the European map inset
# to limit size of the map.Swedish coordinates are shifted slightly south to prevent overlap with  map edge
MetaU <- MetaU %>%
  mutate(PlotRegio = case_when(
      CountryClass == "Africa"       ~ "Africa",
      CountryClass == "Middle East"  ~ "Middle East",
      TRUE                           ~ collection_countryname)) %>%
  mutate(collection_countrylat = case_when(
    CountryClass == "Africa"       ~ 40.5,
    CountryClass == "Middle East"  ~ 40.5,
    collection_countryname == "Sweden" ~ 57.4,
    TRUE                           ~ collection_countrylat))  %>%
  mutate(collection_countrylong = case_when(
    CountryClass == "Africa"       ~ 5.2,
    CountryClass == "Middle East"  ~ 6.8,
    collection_countryname == "Sweden" ~ 14.5,
    TRUE                           ~ collection_countrylong)) 

# Tip labels from each MCC tree
tips_E   <- unique(MCCglE@phylo$tip.label)
tips_NS5 <- unique(MCCglNS5@phylo$tip.label)

# Aggregate sequence counts per country or proxy region for bubble map
MetaU_E <- MetaU %>%
  dplyr::filter(ID %in% tips_E) %>%
  group_by(PlotRegio) %>%
  summarise(n=n(), collection_countrylat = first(collection_countrylat), collection_countrylong = first(collection_countrylong), collection_countryname= first(collection_countryname), CountryClass = first(CountryClass))

MetaU_NS5 <- MetaU %>%
  dplyr::filter(ID %in% tips_NS5) %>%
  group_by(PlotRegio) %>%
  summarise(n=n(), collection_countrylat = first(collection_countrylat), collection_countrylong = first(collection_countrylong), collection_countryname= first(collection_countryname), CountryClass = first(CountryClass))

# Bubble size is scaled proportionally to sequence count
# using the same range across both panels to keep the two maps directly comparable
pts_E <- st_as_sf(MetaU_E,
                coords = c("collection_countrylong","collection_countrylat"), crs = 4326)
pts_NS5 <- st_as_sf(MetaU_NS5,
                  coords = c("collection_countrylong","collection_countrylat"), crs = 4326)

# Create scale jointly for both panels
rng <- range(c(pts_E$n, pts_NS5$n))

# bubble size scaling (area-proportional) 
scale_bubbles <- function(x, from = rng, min_cex = 0.4, max_cex = 4.4) {
  r <- range(x)
  s <- sqrt((x - r[1]) / (r[2] - r[1]))        # area ∝ n
  min_cex + s * (max_cex - min_cex)
}

cex_E   <- scale_bubbles(pts_E$n)
cex_NS5 <- scale_bubbles(pts_NS5$n)

# Set color properties
palette_classes <- names(colors)  # c("Netherlands","Germany","Italy","Europe other","Africa","Middle East")
MetaU$CountryClass <- factor(MetaU$CountryClass, levels = palette_classes)

fill_alpha <- 0.7
fill_cols  <- function(pts) {
  raw_cols <- colors[as.character(pts$CountryClass)]          
  adjustcolor(unname(raw_cols), alpha.f = fill_alpha)
}
border_col <- adjustcolor("#737373", 0.5)

# Function to draw maps
lv <- c(1, 100, 200, 298)
draw_maps <- function(pts, cex_vals, title = NULL, add_legend = FALSE) {
  eval(EUMap)
  plot(
    st_geometry(pts),
    pch = 21,
    cex = cex_vals,
    bg = fill_cols(pts),
    col = border_col,
    lwd = 0.4,
    add = TRUE)
  
  # Add labels for countries with > 20 sequences
  coords <- st_coordinates(pts)
  label_idx <- which(pts$n > 20)
  if (length(label_idx) > 0) {
    text(
      x = coords[label_idx, 1],
      y = coords[label_idx, 2],
      labels = pts$n[label_idx],
      pos = 3,
      cex = 0.7,
      col = colors[as.character(pts$CountryClass[label_idx])],
      font = 2)
    # Function to draw maps
    lv <- c(1, 100, 200, 298)
    draw_maps <- function(pts,cex_vals, title = NULL, add_legend = FALSE) {
      eval(EUMap)
      plot(
        st_geometry(pts),
        pch = 21,
        cex = cex_vals,
        bg = fill_cols(pts),
        col = border_col,
        lwd = 0.4,
        add = TRUE)
      
      # Add labels for countries with > 20 sequences
      coords <- st_coordinates(pts)
      label_idx <- which(pts$n > 30)
      if (length(label_idx) > 0) {
        text(
          x = coords[label_idx, 1],
          y = coords[label_idx, 2],
          labels = pts$n[label_idx],
          pos = 3,
          cex = 0.7,
          col = colors[as.character(pts$CountryClass[label_idx])],
          font = 2)
      }
      
      if (add_legend) {
        legend(
          "topright",
          title = "N USUV sequences",
          legend = lv,
          pch = 21,
          pt.bg  = adjustcolor("grey85", 0.9),
          pt.cex = scale_bubbles(lv),
          col    = border_col,
          pt.lwd = 0.4,
          lty = 0,
          lwd = NA,
          bty = "n")
        usr <- par("usr")
        scalebar(
          d = 500,
          # length of scale bar in km
          xy = c(usr[1] + 0.65 * (usr[2] - usr[1]), usr[3] + 0.12 * (usr[4]-usr[3])),
          type = "bar",
          divs = 2,
          # number of subdivisions
          below = "km",
          lonlat = TRUE,
          cex = 0.7,
          col = "grey85")
      }
    }
  }
}
# Open two panels and draw maps
pdf(paste0(FigRepo,"GlobalDTA/","Maps",".pdf"), 
    width = 8.27, height = 4.13)
op <- par(no.readonly = TRUE); on.exit(par(op))
par(mfrow = c(1,2), mar = c(0,0,0,0), xaxs = "i", yaxs = "i")
draw_maps(pts_E,   cex_E,   title = "nfg+E",add_legend = FALSE)
draw_maps(pts_NS5, cex_NS5, title = "nfg+NS5", add_legend = TRUE)  
dev.off()


## Fig. 2: Time-scaled phylogenies and Effective pop size------------------------
# A:Time-scaled phylogenies

DirSkygrid <- "2_EffectivePopSizeReconstruction/Skygrid_Joint_Clade1_2"

MCCSky1 <- read.beast(file.path(DirSkygrid, "NLCladesJoint_Skygrid.NLClade1_MCC.tree"))
MCCSky2 <- read.beast(file.path(DirSkygrid, "NLCladesJoint_Skygrid.NLClade2_MCC.tree"))

# choose a common y-extent based on the larger tree
y_max <- max(Ntip(MCCSky1), Ntip(MCCSky2)) + 1 

s1 <- ggtree(MCCSky1,mrsd = "2023-03-30", ladderize = TRUE, right = TRUE, size = 0.32, color = "#cc1337") +
  theme_tree2() +
  geom_rootedge(rootedge = 0.8, color = "#cc1337")

s1.1 <- s1  +
  geom_tippoint(color = "#cc1337", size = 1.2)+
  geom_range("height_0.95_HPD", color="#cc1337", size=0.9, alpha=0.15) + 
  geom_nodepoint(aes(subset = as.numeric(posterior) > 0.9), 
                 shape = 16, size = 1, color = "#cc1337") +  
  Arbo_theme_grid() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),      
        axis.ticks.x = element_blank(), 
        axis.text.y = element_blank(),      
        axis.ticks.y = element_blank(),     
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank()) +
  ggtitle(NULL, subtitle = "Clade 1") +
  scale_x_continuous(breaks = seq(2010, 2024, by = 2),limits = c(2009.5, 2024)) +
  coord_cartesian(ylim = c(-1.7, y_max), expand = FALSE)

s2 <- ggtree(MCCSky2,mrsd = "2023-09-28", ladderize = TRUE, right = TRUE, size = 0.32, color = "#cc1337") +
  theme_tree2() +
  geom_rootedge(rootedge = 0.8, color = "#cc1337")

s2.1 <- s2  +
  geom_tippoint(color = "#cc1337", size = 1.2)+
  geom_range("height_0.95_HPD", color="#cc1337", size=0.9, alpha=0.15) + 
  geom_nodepoint(aes(subset = as.numeric(posterior) > 0.9), 
                 shape = 16, size = 1, color = "#cc1337") +  
  Arbo_theme_grid() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8.89, color = "#231f20", angle = 45, hjust = 1), # X-axis text color
        axis.text.y = element_blank(),      
        axis.ticks.y = element_blank(),     
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank()) +
  ggtitle(NULL, subtitle = "Clade 2") +   # Subtitle only 
  scale_x_continuous(breaks = seq(2010, 2024, by = 2),limits = c(2009.5, 2024)) +
  coord_cartesian(ylim = c(-1.7, y_max), expand = FALSE)


s <- s1.1/s2.1

# B:Ne
# Define grid points as they were set up for the skygrid analysis
# Meteorological seasons
# winter: december, januari, februari
# spring: maart, april, mei
# summer: juni, juli, augustus
# autumn: september, oktober, november
# points are placed at the start of
# each meteorological season.
grid.points.spring <- paste0(2010:2023, c("-03-01"))
grid.points.summer <- paste0(2010:2023, c("-06-01"))
grid.points.fall <- paste0(2010:2023, c("-09-01"))
grid.points.winter <- paste0(2010:2023, c("-12-01"))

#merge:
separator <- "XXXX"

combo <- paste(grid.points.spring, grid.points.summer, grid.points.fall,grid.points.winter,
               sep = separator)
grid.points <- unlist(strsplit(x = combo, split = separator))

# Truncate grid at the most recent sampling date of each clade
grid.points1 <- c(head(grid.points, -3)) #In clade 1 most recent sample in spring -> remove Summer, Autumn, Winter 2023 as we'll use the MRDS
grid.points2 <- c(head(grid.points, -1)) #In clade 2 most recent sample in autumn

grid.points1.decimal <- lubridate::decimal_date(as.Date(grid.points1))
grid.points2.decimal <- lubridate::decimal_date(as.Date(grid.points2))

MRSD1 <- lubridate::decimal_date(as.Date("2023-03-30")) #date of the most recent sample Clade 1: 2023.241
MRSD2 <- lubridate::decimal_date(as.Date("2023-09-28")) #date of the most recent sample Clade 2: 2023.74

# Append MRSD and reverse to run from present to past, matching BEAST output ordering
grid.points1.decimalMRSD1 <- c(grid.points1.decimal, MRSD1)
grid.points1.decimalMRSD1<- rev(grid.points1.decimalMRSD1)

# Ne is estimated between grid points; midpoints of each interval are used as
# the time coordinate for plotting
midpoints1 <- (head(grid.points1.decimalMRSD1, -1) + tail(grid.points1.decimalMRSD1, -1)) / 2
midpoints1 <- c(midpoints1, NA) #pad with NA for the last element

# Same for clade 2
grid.points2.decimalMRSD2 <- c(grid.points2.decimal, MRSD2)
grid.points2.decimalMRSD2<- rev(grid.points2.decimalMRSD2)

# Ne is estimated between grid points; midpoints of each interval are used as
# the time coordinate for plotting
midpoints2 <- (head(grid.points2.decimalMRSD2, -1) + tail(grid.points2.decimalMRSD2, -1)) / 2
midpoints2 <- c(midpoints2, NA) #pad with NA for the last element

# Directory containing BEAST skygrid log files and BNPR-PS output from seraphim
DirectorySky = "2_EffectivePopSizeReconstruction/Skygrid_Joint_Clade1_2/"
DirectoryPref = "2_EffectivePopSizeReconstruction/BNPR_PS_Joint_Clade1_2/Results/"

# Read skygrid logPopSize parameters
logSky <- read_delim(paste0(DirectorySky,"NLCladesJoint_Skygrid.log"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 4)
burnin_fraction <- 0.1
burnInSky <- round(burnin_fraction * dim(logSky)[1]) + 1

# keep only pop.size values
which(colnames(logSky)=="NLClade1.skygrid.logPopSize1") #13
which(colnames(logSky)=="NLClade2.skygrid.logPopSize56") #122

logSky <- logSky[,c(1,13:122)]
logSky1 = logSky[burnInSky:dim(logSky)[1],]
selected_colsSky <- grep("skygrid.logPopSize", colnames(logSky1), value = TRUE)

# Compute posterior mean and 95% HPD interval for each Ne(t) interval;
# Initialize empty vectors to store results
names_vector <- c()
mean_vector <- c()
low_hpd_vector <- c()
high_hpd_vector <- c()

# Loop through each selected column: mean value, hpd interval -> 
for (col_name in selected_colsSky) {
  mean_val <- mean(logSky1[[col_name]])
  hpd_low <- HDInterval::hdi(logSky1[[col_name]])[1]
  hpd_high <-  HDInterval::hdi(logSky1[[col_name]])[2]
  #store results in the vectors
  names_vector <- c(names_vector, col_name)
  mean_vector <- c(mean_vector, mean_val)
  low_hpd_vector <- c(low_hpd_vector, hpd_low)
  high_hpd_vector <- c(high_hpd_vector,hpd_high)
}

# Create a data frame from the collected results
result_dfSky <- data.frame(
  name = names_vector,
  mean = mean_vector,
  low_HPD = low_hpd_vector,
  high_HPD = high_hpd_vector
)

# Read BNPR-PS log
logPref <- read_delim(paste0(DirectoryPref,"USUV.NL.Mainclade_PrefSampling_Seasons_6.pop"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 3)
str(logPref)
burnInPref <- round(burnin_fraction * dim(logPref)[1]) +1

#keep only pop.size values
which(colnames(logPref)=="NLClade1.logPop1") #2
which(colnames(logPref)=="NLClade1.logPop54") #55
which(colnames(logPref)=="NLClade2.logPop1") #110
which(colnames(logPref)=="NLClade2.logPop56") #165

logPref <- logPref[,c(1:55, 110:165)] #keep only pop.size values
logPref1 = logPref[burnInPref:dim(logPref)[1],]
selected_colsPref <- grep("logPop", colnames(logPref1), value = TRUE)

# Compute posterior mean and 95% HPD interval for each Ne(t) interval;
# Initialize empty vectors to store results
names_vector <- c()
mean_vector <- c()
low_hpd_vector <- c()
high_hpd_vector <- c()

for (col_name in selected_colsPref) {
  mean_val <- mean(logPref1[[col_name]])
  hpd_low <- HDInterval::hdi(logPref1[[col_name]])[1]
  hpd_high <-  HDInterval::hdi(logPref1[[col_name]])[2]
  #store results
  names_vector <- c(names_vector, col_name)
  mean_vector <- c(mean_vector, mean_val)
  low_hpd_vector <- c(low_hpd_vector, hpd_low)
  high_hpd_vector <- c(high_hpd_vector,hpd_high)
}

result_dfPref <- data.frame(
  name = names_vector,
  mean = mean_vector,
  low_HPD = low_hpd_vector,
  high_HPD = high_hpd_vector
)

result_dfPref$method <- "Pref"
result_dfSky$method <- "Sky"

# Bind posterior Ne estimates to the midpoint of each seasonal interval,
# then bind skygrid and BNPR-PS results
result_dfSky1 <- result_dfSky %>%
  filter(grepl("NLClade1", name)) %>%
  cbind(midpoints1)

result_dfPref1 <- result_dfPref %>%
  filter(grepl("NLClade1", name)) %>%
  cbind(midpoints1)

dfCombined1 <- rbind(result_dfPref1, result_dfSky1)
colnames(dfCombined1) <- c("name","mean","low_HPD","high_HPD","method","time")

result_dfSky2 <- result_dfSky %>%
  filter(grepl("NLClade2", name)) %>%
  cbind(midpoints2)

result_dfPref2 <- result_dfPref %>%
  filter(grepl("NLClade2", name)) %>%
  cbind(midpoints2)

dfCombined2 <- rbind(result_dfPref2, result_dfSky2)
colnames(dfCombined2) <- c("name","mean","low_HPD","high_HPD","method","time")

# Yearly vertical lines are positioned at the start of spring (decimal year ≈ .16)
years <- grid.points2.decimalMRSD2[
  grepl("\\.16$", sprintf("%.2f", grid.points2.decimalMRSD2))
]

# Define breaks
breaks <- grid.points2.decimalMRSD2

# Plot Ne(t) time series comparing skygrid (blue) and BNPR-PS (red) for Clade 1
E1 <- dfCombined1 %>%
  ggplot(aes(x = time, y = mean, color = method, fill = method)) +
  
  # Add reference lines
  geom_hline(yintercept = 0, color = "#231f20") +  
  geom_vline(xintercept = years, color = "#d1d3d4", linewidth = 0.1) +  # Vertical lines marking years (start of spring)
  
  # Sky method: credible interval (ribbon) and mean line
  geom_ribbon(
    data = subset(dfCombined1, method == "Sky"),
    aes(ymin = low_HPD, ymax = high_HPD, fill = method),
    alpha = 0.2, 
    color = NA    
  ) +
  geom_line(
    data = subset(dfCombined1, method == "Sky"),
    aes(color = method)
  ) +
  
  # Pref method: credible interval (ribbon) and mean line
  geom_ribbon(
    data = subset(dfCombined1, method == "Pref"),
    aes(ymin = low_HPD, ymax = high_HPD, fill = method),
    alpha = 0.2,  
    color = NA) +
  geom_line(
    data = subset(dfCombined1, method == "Pref"),
    aes(color = method)) +
  
  # Set custom colors for each method
  scale_color_manual(values = c("Sky" = "#006eb8", "Pref" = "#cc1337")) +
  scale_fill_manual(values = c("Sky" = "#006eb8", "Pref" = "#cc1337")) +
  
  scale_x_continuous(
    limits = c(2011.162, 2023.703),
    expand = c(0, 0),  
    breaks = years,
    labels = NULL) +
  coord_cartesian(
    xlim = c(2011.414, 2023.703),  
    expand = FALSE) +
    scale_y_continuous(
    limits = c(-9.5, 13.75),
    expand = c(0, 0)) +
  
  labs(
    x = NULL,  
    y = expression(log(N[e]))) +
  
  Arbo_theme_axes() +
  
  # Configure legend
  guides(
    fill = guide_legend(title = NULL),   
    color = guide_legend(title = NULL)) +
  theme(
    axis.ticks.x = element_blank(), 
    legend.position = c(0.8, 0.001),     
    legend.justification = c(0, 0),       
    legend.background = element_blank(),  
    legend.key = element_blank())

# Plot Ne(t) time series comparing skygrid (blue) and BNPR-PS (red) for Clade 2
E2 <- dfCombined2 %>%
  ggplot(aes(x = time, y = mean, color = method, fill = method)) +
  
  # Add reference lines
  geom_hline(yintercept = 0, color = "#231f20") +  
  geom_vline(xintercept = years, color = "#d1d3d4", linewidth = 0.1) +  # Vertical lines marking years (start of spring)
  
  # Sky method: credible interval (ribbon) and mean line
  geom_ribbon(
    data = subset(dfCombined2, method == "Sky"),
    aes(ymin = low_HPD, ymax = high_HPD, fill = method),
    alpha = 0.2,
    color = NA) +
  geom_line(
    data = subset(dfCombined2, method == "Sky"),
    aes(color = method)) +
  
  # Pref method: credible interval (ribbon) and mean line
  geom_ribbon(
    data = subset(dfCombined2, method == "Pref"),
    aes(ymin = low_HPD, ymax = high_HPD, fill = method),
    alpha = 0.2,
    color = NA) +
  geom_line(
    data = subset(dfCombined2, method == "Pref"),
    aes(color = method)) +
  
  # Set custom colors for each method
  scale_color_manual(values = c("Sky" = "#006eb8", "Pref" = "#cc1337")) +
  scale_fill_manual(values = c("Sky" = "#006eb8", "Pref" = "#cc1337")) +
  
  scale_x_continuous(
    limits = c(2011.162, 2023.703),
    expand = c(0, 0),  
    breaks = years,
    labels = as.character(floor(years))) +
  coord_cartesian(
    xlim = c(2011.414, 2023.703),  # Zoom to specific x-range
    expand = FALSE) +
  
  scale_y_continuous(
    limits = c(-9.5, 13.75),
    expand = c(0, 0)) +
  
  labs(
    x = NULL,  # No x-axis label
    y = expression(log(N[e]))) +
  
  Arbo_theme_axes() +
  
  guides(
    fill = guide_legend(title = NULL),
    color = guide_legend(title = NULL)) +
  theme(
    legend.position = c(0.8, 0.001),   
    legend.justification = c(0, 0),      
    legend.background = element_blank(), 
    legend.key = element_blank()
  )

E1 <- E1 + ggtitle(NULL, subtitle = "Clade 1")
E2 <- E2 + ggtitle(NULL, subtitle = "Clade 2")  

E1/E2

# Heatmap overlay: seasonal sequencing counts (per clade) and USUV-positive dead bird
# detections are displayed as a two-row heatmap below each Ne(t) panel to show the
# temporal sampling covariate used in the BNPR-PS model

# Load count data
df_counts <- read.csv(paste0("2_EffectivePopSizeReconstruction/BNPR_PS_Joint_Clade1_2/","NUSUVDeadBirds_SeqMainNLClades_Season.csv"),
                      header = TRUE,
                      sep = ",",
                      stringsAsFactors = FALSE)

# Reverse row order to align seasonal counts with the time axis used in the Ne plots
rev_counts <- df_counts[nrow(df_counts):1, ]
rev_counts$midpoints <- midpoints2

# Assign a default midpoint (2010.04) for intervals before first skygrid point
rev_counts <- rev_counts %>%
  mutate(midpoints = ifelse(is.na(midpoints), 2010.04, midpoints))

# Select relevant columns for clade 1
CountsC1 <- rev_counts %>%
  dplyr::select(midpoints, NSeqs_NLClade1, NUSUVDeadBirds)

# Reshape to long format for tile-based heatmap rendering
CountsC1_long <- CountsC1 %>%
  pivot_longer(
    cols = -midpoints,
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(
    midpoints = as.numeric(midpoints),
    value = as.numeric(value),
    variable = factor(variable)
  )

# Get number of variables (rows in heatmap)
n1 <- nlevels(CountsC1_long$variable)

# Create color scale: white-to-red gradient;
# the colour bar range is shared across both Clade panels to allow direct comparison

rose_pal <- c("#f8e9e4", "#cc1337")
col <- c("white", colorRampPalette(rose_pal)(99))

# Define heatmap position 
band_top    <- -7.5
band_bottom <- -9.5

# Calculate tile positions for each variable
CountsC1_long2 <- CountsC1_long %>%
  mutate(
    row = case_when(
      variable == "NSeqs_NLClade1"   ~ 2,
      variable == "NUSUVDeadBirds"   ~ 1),
    tile_y = band_bottom + (row - 0.5) * (band_top - band_bottom) / n1,
    tile_h = (band_top - band_bottom) / n1)

# Add heatmap overlay to clade 1 Ne(t) panel
E1_overlay <- E1 +
  # Enable new fill scale for heatmap (independent of main plot)
  ggnewscale::new_scale_fill() +
  
  # Add heatmap tiles
  geom_tile(
    data = CountsC1_long2,
    aes(x = midpoints, y = tile_y, fill = value),
    inherit.aes = FALSE,
    width = 0.25,        # Tile width in time units
    height = 1,        # Tile height on y-axis
    color = "#e5e5e5",   # Light gray border
    linewidth = 0.15
  ) +
  
  # Configure heatmap color gradient
  scale_fill_gradientn(
    colours = col,
    limits = c(0, max(CountsC1_long2$value, na.rm = TRUE)),
    oob = scales::squish,  # Handle values outside limits
    guide = "none"         # Hide color scale legend
  ) +
  
  # Frame the heatmap with reference lines
  geom_hline(yintercept = band_top, color = "#231f20", linewidth = 0.3) +
  geom_hline(yintercept = band_top - 1, color = "#737373", linewidth = 0.1) +
  
  theme(
    legend.position = c(0.8, 0.001),
    legend.justification = c(0, 0),
    legend.background = element_blank(),
    legend.key = element_blank()
  ) +
  
  # Add row labels for heatmap
  annotate(
    "text",
    x = 2023, y = band_top - 0.5,
    label = "Sequenced (Clade 1)",
    hjust = 0, size = 3,
  ) +
  annotate(
    "text",
    x = 2023, y = band_bottom + 0.5,
    label = "USUV+ birds",
    hjust = 0, size = 3,
  ) 

# Select relevant columns for clade 2
CountsC2 <- rev_counts %>%
  dplyr::select(midpoints, NSeqs_NLClade2, NUSUVDeadBirds)

# Convert to long format for plotting
CountsC2_long <- CountsC2 %>%
  pivot_longer(
    cols = -midpoints,
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(
    midpoints = as.numeric(midpoints),
    value = as.numeric(value),
    variable = factor(variable)
  )

# Get number of variables (rows in heatmap)
n2 <- nlevels(CountsC2_long$variable)

# Create color scale: white-to-red gradient;
# the colour bar range is shared across both Clade panels to allow direct comparison
rose_pal <- c("#f8e9e4", "#cc1337")
col <- c("white", colorRampPalette(rose_pal)(99))

# Define heatmap position on y-axis
band_top    <- -7.5  # Top edge of heatmap
band_bottom <- -9.5   # Bottom edge of heatmap

# Calculate tile positions for each variable
CountsC2_long2 <- CountsC2_long %>%
  mutate(
    row = case_when(
      variable == "NSeqs_NLClade2"   ~ 2,
      variable == "NUSUVDeadBirds"   ~ 1),
    tile_y = band_bottom + (row - 0.5) * (band_top - band_bottom) / n1,
    tile_h = (band_top - band_bottom) / n1
  )

# Add heatmap overlay to clade 2 Ne(t) panel
E2_overlay <- E2 +
  # Enable new fill scale for heatmap (independent of main plot)
  ggnewscale::new_scale_fill() +
  
  # Add heatmap tiles
  geom_tile(
    data = CountsC2_long2,
    aes(x = midpoints, y = tile_y, fill = value),
    inherit.aes = FALSE,
    width = 0.25,        # Tile width in time units
    height = 1,        # Tile height on y-axis
    color = "#e5e5e5",   # Light gray border
    linewidth = 0.15
  ) +
  
  # Configure heatmap color gradient
  scale_fill_gradientn(
    colours = col,
    limits = c(0, max(CountsC1_long2$value, na.rm = TRUE)),
    oob = scales::squish,  # Handle values outside limits
    guide = guide_colorbar(title = "Count", title.position = "top")  # Show continuous legend
  ) +
  
  # Frame the heatmap with reference lines
  geom_hline(yintercept = band_top, color = "#231f20", linewidth = 0.3) +
  geom_hline(yintercept = band_top - 1, color = "#737373", linewidth = 0.1) +
  
  # Maintain legend settings from original plot
  guides(
    color = guide_legend(title = NULL)
  ) +
  theme(
    legend.position = c(0.8, 0.001),
    legend.justification = c(0, 0),
    legend.background = element_blank(),
    legend.key = element_blank()
  ) +
  
  # Add row labels for heatmap
  annotate(
    "text",
    x = 2023, y = band_top - 0.5,
    label = "Sequenced (Clade 2)",
    hjust = 0, size = 3
  ) +
  annotate(
    "text",
    x = 2023, y = band_bottom + 0.5,
    label = "USUV+ birds",
    hjust = 0, size = 3
  )

# Season abbreviations are overlaid on the heatmap band for a single year (2012)
# Define label positions 
season_labs <- data.frame(
  season = c("SP", "SU", "AU","WI"),  # Winter, Spring, Summer, Autumn
  x = c(2012.290, 2012.541, 2012.791,2013.038)
)

# Create final panel with season annotations
E2_overlay.2 <- E2_overlay +
  annotate(
    "text",
    x = season_labs$x,
    y = -6.5,
    label = season_labs$season,
    angle = 45,
    vjust = 0.5,
    size = 2,
    color = "grey35"
  )

Ex <- E1_overlay/E2_overlay.2

pdf(paste0(FigRepo,"TwoCladesTimeScaled/","SkygridBNPRPS",".pdf"), 
    width = 8.15, height = 5.71)  # A5 with small offset
s|Ex
dev.off()

# Export portrait A5
pdf(paste0(FigRepo,"TwoCladesTimeScaled/","Trees",".pdf"), 
    width = 4.13, height = 5.83)
s1.1/s2.1
dev.off()

## Fig. 3: Continuous phylogeography----------------------------------------

DirCPhyloGeo <- "3_ContinuousPhyloGeo/"
tree.files <- list.files(DirCPhyloGeo, pattern = ".trees")

mcc_tab1 <- read.csv(file.path(DirCPhyloGeo, "250812_NLClade1_SkyRRW.csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)
mcc_tab2 <- read.csv(file.path(DirCPhyloGeo, "250812_NLClade2_SkyRRW.csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)
mcc_tabj <- read.csv(file.path(DirCPhyloGeo, "All_clades.csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)

minYearv = min(mcc_tab2[,"startYear"])

localTreesDirectory1 <- file.path(DirCPhyloGeo, "250812_NLClade1_SkyRRW_ext")
localTreesDirectory2 <- file.path(DirCPhyloGeo, "250812_NLClade2_SkyRRW_ext")  

# Estimate the HPD region for each time slice
startDatum1 <- suppressWarnings(min(as.numeric(mcc_tab1$startYear), na.rm = TRUE))
startDatum2 <- suppressWarnings(min(as.numeric(mcc_tab2$startYear), na.rm = TRUE))

percentage <- 80
prob <- percentage / 100
precision <- 1          # years; use 1/12 for monthly slices
nberOfExtractionFiles <- 900  # ensure this matches your available trees

# spreadGraphic2 produces a list of distinct spatial polygon data frames, with one data frame for each time slice
# the estiamtion of the HPD region is based on all the ending positions of phylogenetic branches whose ending time falls within the considered time slice.
polygons1 <- spreadGraphic2(localTreesDirectory1, nberOfExtractionFiles, prob, 2011, precision)
polygons2 <- spreadGraphic2(localTreesDirectory2, nberOfExtractionFiles, prob, 2010, precision)

# Visualize Spread - Co-plot the HPD regions and MCC tree, per clade

# Read shapefile for background map (Provinces of the Netherlands in WGS 1984)
admin1 <- shapefile("Shapefiles/NL_ADM0_Land.shp")

# Set up colour palette
full <- colorRampPalette(c("#00375c", "#006eb8","#d5abbd", "#cc1337", "#660a1c"))(141)

blue_part <- full[25:67]   
red_part  <- full[73:130]  

# Combine (blue on left, red on right)
colourScale <- c(blue_part, red_part)
length(colourScale) #For the following script, make sure colour scale lenght is 101

barplot(rep(1, length(colourScale)), col = colourScale, border = NA, space = 0,
        main = "Discrete Color Scale Visualization",
        xlab = "Color Index", ylab = "", axes = FALSE)

# Dummy raster (just to give a scale range)
rast <- raster(matrix(nrow=1, ncol=2))
rast[1] <- minYear
rast[2] <- maxYear

# Plot the legend only
plot(rast,
  legend.only   = TRUE,
  col           = colourScale,horizontal    = TRUE, legend.width  = 0.5,legend.shrink = 0.3,
  smallplot     = c(0.15, 0.85, 0.1, 0.2),
  legend.args   = list(text = "", cex = 0.7, line = 0.3, col = "gray30"),
  axis.args     = list(cex.axis = 0.75,lwd      = 0,lwd.tick = 0.2,col.tick = "gray30",tck      = -0.9,col      = "gray30",col.axis = "gray30",line     = 0,mgp      = c(0, 0.1, 0)
  )
)

minYear = 2010
maxYear = 2024

# Map each polygon's time coordinate to a colour index so that the HPD polygon
# and MCC branch/node colours share the same temporal colour scale
rast <- raster(matrix(nrow=1, ncol=length(colourScale)))
values(rast) <- seq(minYear, maxYear, length.out=length(colourScale))

polygons_colours1 = rep(NA, length(polygons1))
for (i in 1:length(polygons1)){
  date = as.numeric(names(polygons1[[i]]))
  polygon_index1 = round((((date-minYear)/(maxYear-minYear))*100)+1)
  polygons_colours1[i] = paste0(colourScale[polygon_index1],"50") #The value "40" represents a level of transparency
}

polygons_colours2 = rep(NA, length(polygons2))
for (i in 1:length(polygons2)){
  date = as.numeric(names(polygons2[[i]]))
  polygon_index2 = round((((date-minYear)/(maxYear-minYear))*100)+1)
  polygons_colours2[i] = paste0(colourScale[polygon_index2],"50") #The value "40" represents a level of transparency
}

# Spread visualizations - Eight-panel layout: 4 time periods × 2 clades. 
# The HPD regions and MCC tree are co-plotting, per clade; wihtin time period
pdf(paste0(DirCPhyloGeo,"RRW_dispersal_plot/","RRW_temporalspread_test","_v2.pdf"), 
    width = 8.3, height = 5.8)

par(mfrow=c(2,4), mar=c(0.5,0.5,0.5,0.5), oma=c(1,1,1,1), mgp=c(0,0.4,0), lwd=0.1, bty="o")

#Time cut-offs for panels
cutOffs = c(2017, 2019, 2022, 2024); croppingPolygons = FALSE

# ---------- Clade 1
for (h in 1:length(cutOffs)){
  #Plot administrative boundaries polygons
  par(mar = c(0.5,0.5,0.5,0.5), xaxs = "i", yaxs = "i")
  eval(BaseMap.Riv)
  #Dummy raster to drive the legend scale from startDatum to max endYear
  rast = raster(matrix(nrow=1, ncol=2)); rast[1] = minYear; rast[2] = maxYear
  
  #Plot HPD polygons in the intervals between cutoffs
  for (i in 1:length(polygons1)){
    val <- as.numeric(names(polygons1[[i]]))
    #everything before the first cutoff
    if (h == 1) {       
      if (val < cutOffs[1]) {
        pol = polygons1[[i]]; crs(pol) = crs(admin1)
        if (croppingPolygons == TRUE) pol <- crop(pol, admin1)
        plot(pol, axes=F, col=polygons_colours1[[i]], add=TRUE, border=NA)
        }
      } else {
      #everything between cutoff[h-1] and cutoff[h]
      if (val < cutOffs[h] & val >= cutOffs[h-1]) 
      {
        pol = polygons1[[i]]; crs(pol) = crs(admin1)
        if (croppingPolygons == TRUE) pol <- crop(pol, admin1)
        plot(pol, axes=F, col=polygons_colours1[[i]], add=TRUE, border=NA)
      }
      }
    }
  
  #Plot MCC branches whose endYear falls inside the current cutoff interval
  for (i in 1:nrow(mcc_tab1)) {
    if (h == 1) {
      if (mcc_tab1[i,"endYear"] < cutOffs[h])
        {
      curvedarrow(cbind(mcc_tab1[i,"startLon"],mcc_tab1[i,"startLat"]), cbind(mcc_tab1[i,"endLon"],mcc_tab1[i,"endLat"]), arr.length=0,
                  arr.width=0, lwd=0.1, lty=1, lcol="gray30", arr.col=NA, arr.pos=F, curve=0.1, dr=NA, endhead=F)
      }
    } else {
      #everything between cutoff[h-1] and cutoff[h]
      if (mcc_tab1[i,"endYear"]  < cutOffs[h] & mcc_tab1[i,"endYear"] >= cutOffs[h-1]) 
      {
        curvedarrow(cbind(mcc_tab1[i,"startLon"],mcc_tab1[i,"startLat"]), cbind(mcc_tab1[i,"endLon"],mcc_tab1[i,"endLat"]), arr.length=0,
                    arr.width=0, lwd=0.1, lty=1, lcol="gray30", arr.col=NA, arr.pos=F, curve=0.1, dr=NA, endhead=F)
      }
    }
  }
  #Plot MCC nodes (points) whose endYear falls inside the current cutoff interval
  for (i in 1:nrow(mcc_tab1)) {
    if (h == 1) {
      if (mcc_tab1[i,"endYear"] < cutOffs[h])
      {
        startYears_index = (((mcc_tab1[i,"startYear"]-minYear)/(maxYear-minYear))*100)+1
        points(mcc_tab1[i,"startLon"], mcc_tab1[i,"startLat"], pch=16, col=colourScale[startYears_index], cex=0.7)
        points(mcc_tab1[i,"startLon"], mcc_tab1[i,"startLat"], pch=1, col="gray30", lwd=0.1, cex=0.7)
      }
    } else {
      if (mcc_tab1[i,"endYear"] < cutOffs[h] & mcc_tab1[i,"endYear"] >= cutOffs[h-1])
      {
        startYear_index = (((mcc_tab1[i,"startYear"]-minYear)/(maxYear-minYear))*100)+1
        points(mcc_tab1[i,"startLon"], mcc_tab1[i,"startLat"], pch=16, col=colourScale[startYear_index], cex=0.7)
        points(mcc_tab1[i,"startLon"], mcc_tab1[i,"startLat"], pch=1, col="gray30", lwd=0.1, cex=0.7)
      }
    }
  }
  for (i in 1:nrow(mcc_tab1)) {
    if (h == 1) {
      if (mcc_tab1[i,"endYear"] < cutOffs[h])
      {
        endYear_index = (((mcc_tab1[i,"endYear"]-minYear)/(maxYear-minYear))*100)+1
        points(mcc_tab1[i,"endLon"], mcc_tab1[i,"endLat"], pch=16, col=colourScale[endYear_index], cex=0.7)
        points(mcc_tab1[i,"endLon"], mcc_tab1[i,"endLat"], pch=1, col="gray30", lwd=0.2, cex=0.7)
      }
    } else {
      if (mcc_tab1[i,"endYear"]  < cutOffs[h] & mcc_tab1[i,"endYear"] >= cutOffs[h-1]) 
      {
        endYear_index = (((mcc_tab1[i,"endYear"]-minYear)/(maxYear-minYear))*100)+1
        points(mcc_tab1[i,"endLon"], mcc_tab1[i,"endLat"], pch=16, col=colourScale[endYear_index], cex=0.7)
        points(mcc_tab1[i,"endLon"], mcc_tab1[i,"endLat"], pch=1, col="gray30", lwd=0.1, cex=0.7)
      }
    }
  }
  
  if (h == 1) {
    label <- paste0(floor(startDatum1), "-", cutOffs[h] - 1)
  } else {
    label <- paste0(cutOffs[h - 1], "-", cutOffs[h] - 1)
  }
  usr <- par("usr")  #returns c(xmin, xmax, ymin, ymax)
  text(x = usr[1] + 0.05*(usr[2]-usr[1]),  # xmin plus 5% width
       y = usr[4] - 0.1*(usr[4]-usr[3]),  # ymax minus 10% height
       labels =  label,
       adj = c(0,1),
       col="gray30") # upper left
}
# ---------- Clade 2 
for (h in 1:length(cutOffs)){
  eval(BaseMap.Riv)
  
  #Plot HPD polygons in the intervals between cutoffs
  for (i in 1:length(polygons2)){
    val <- as.numeric(names(polygons2[[i]]))
    #everything before the first cutoff
    if (h == 1) {       
      if (val < cutOffs[1]) {
        pol = polygons2[[i]]; crs(pol) = crs(admin1)
        if (croppingPolygons == TRUE) pol <- crop(pol, admin1)
        plot(pol, axes=F, col=polygons_colours2[[i]], add=TRUE, border=NA)
      }
    } else {
      #everything between cutoff[h-1] and cutoff[h]
      if (val < cutOffs[h] & val >= cutOffs[h-1]) 
      {
        pol = polygons2[[i]]; crs(pol) = crs(admin1)
        if (croppingPolygons == TRUE) pol <- crop(pol, admin1)
        plot(pol, axes=F, col=polygons_colours2[[i]], add=TRUE, border=NA)
      }
    }
    }
  
  #Plot MCC branches whose endYear falls inside the current cutoff interval
  for (i in 1:nrow(mcc_tab2)) {
    if (h == 1) {
      if (mcc_tab2[i,"endYear"] < cutOffs[h])
      {
        curvedarrow(cbind(mcc_tab2[i,"startLon"],mcc_tab2[i,"startLat"]), cbind(mcc_tab2[i,"endLon"],mcc_tab2[i,"endLat"]), arr.length=0,
                    arr.width=0, lwd=0.2, lty=1, lcol="gray30", arr.col=NA, arr.pos=F, curve=0.1, dr=NA, endhead=F)
      }
    } else {
      #everything between cutoff[h-1] and cutoff[h]
      if (mcc_tab2[i,"endYear"]  < cutOffs[h] & mcc_tab2[i,"endYear"] >= cutOffs[h-1]) 
      {
        curvedarrow(cbind(mcc_tab2[i,"startLon"],mcc_tab2[i,"startLat"]), cbind(mcc_tab2[i,"endLon"],mcc_tab2[i,"endLat"]), arr.length=0,
                    arr.width=0, lwd=0.2, lty=1, lcol="gray30", arr.col=NA, arr.pos=F, curve=0.1, dr=NA, endhead=F)
      }
    }
  }
  #Plot MCC nodes (points) whose endYear falls inside the current cutoff interval
  for (i in 1:nrow(mcc_tab2)) {
    if (h == 1) {
      if (mcc_tab2[i,"endYear"] < cutOffs[h])
      {
        startYears_index = (((mcc_tab2[i,"startYear"]-minYear)/(maxYear-minYear))*100)+1
        points(mcc_tab2[i,"startLon"], mcc_tab2[i,"startLat"], pch=16, col=colourScale[startYears_index], cex=0.7)
        points(mcc_tab2[i,"startLon"], mcc_tab2[i,"startLat"], pch=1, col="gray30", lwd=0.2, cex=0.7)
      }
    } else {
      if (mcc_tab2[i,"endYear"] < cutOffs[h] & mcc_tab2[i,"endYear"] >= cutOffs[h-1])
      {
        startYears_index = (((mcc_tab2[i,"startYear"]-minYear)/(maxYear-minYear))*100)+1
        points(mcc_tab2[i,"startLon"], mcc_tab2[i,"startLat"], pch=16, col=colourScale[startYears_index], cex=0.7)
        points(mcc_tab2[i,"startLon"], mcc_tab2[i,"startLat"], pch=1, col="gray30", lwd=0.2, cex=0.7)
      }
    }
  }
  for (i in 1:nrow(mcc_tab2)) {
    if (h == 1) {
      if (mcc_tab2[i,"endYear"] < cutOffs[h])
      {
        endYear_index = (((mcc_tab2[i,"endYear"]-minYear)/(maxYear-minYear))*100)+1
        points(mcc_tab2[i,"endLon"], mcc_tab2[i,"endLat"], pch=16, col=colourScale[endYear_index], cex=0.7)
        points(mcc_tab2[i,"endLon"], mcc_tab2[i,"endLat"], pch=1, col="gray30", lwd=0.2, cex=0.7)
      }
    } else {
      if (mcc_tab2[i,"endYear"]  < cutOffs[h] & mcc_tab2[i,"endYear"] >= cutOffs[h-1]) 
      {
        endYear_index = (((mcc_tab2[i,"endYear"]-minYear)/(maxYear-minYear))*100)+1
        points(mcc_tab2[i,"endLon"], mcc_tab2[i,"endLat"], pch=16, col=colourScale[endYear_index], cex=0.7)
        points(mcc_tab2[i,"endLon"], mcc_tab2[i,"endLat"], pch=1, col="gray30", lwd=0.2, cex=0.7)
      }
    }
  }
  
  if (h == 1) {
    label <- paste0(floor(startDatum2), "-", cutOffs[h] - 1)
        } else {
      label <- paste0(cutOffs[h - 1], "-", cutOffs[h] - 1)
      }
  usr <- par("usr")  #returns c(xmin, xmax, ymin, ymax)
  text(x = usr[1] + 0.05*(usr[2]-usr[1]),  # xmin plus 5% width
       y = usr[4] - 0.1*(usr[4]-usr[3]),  # ymax minus 10% height
       labels =  label,
       adj = c(0,1),
       col="gray30") # upper left
  
  #Plot legend: only on the last panel
  if (h == length(cutOffs)){
    # ---- Time legend (left) ----
    plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, 
         smallplot=c(0.510,0.870,0.100,0.115),
         legend.args=list(text="", cex=0.7, line=0.3, col="gray30"), horizontal=T,
         axis.args=list(cex.axis=0.75, lwd=0, lwd.tick=0.2, col.tick="gray30", tck=-0.9, col="gray30", col.axis="gray30", line=0, mgp=c(0,0.10,0)))
    
    # ---- Scalebar (right) ----
    usr <- par("usr")  # c(xmin, xmax, ymin, ymax)
    scalebar(d = 100,        # length of scale bar in km
             xy = c(usr[1] + 0.05*(usr[2]-usr[1]), usr[3] + 0.2*(usr[4]-usr[3])),
             type = "bar", 
             divs = 4,      # number of subdivisions
             below = "km", 
             lonlat = TRUE,
             cex=0.7, col="gray30")
  }
}  

dev.off()

## Fig. 4: Discrete phylogeography ------------------------------------------

# Extract the information from MCC tree and posterior trees (1001 total  trees)
DirProvDTA = "4_DiscretePhyloGeo_Provinces_PerClade/"
DTA_analysis = "DTA_EmpTrees/"
TSW_analysis <- "DTA_EmpTrees_TSW/"

nberOfTreesToSample = 900 
nberOfExtractionFiles <- nberOfTreesToSample
log <- read_delim(paste0(DirProvDTA, DTA_analysis,"Clade1_ProvincesDTA_EmpTrees.log"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 4)
burnin_fraction <- 0.1
burnIn <- round(burnin_fraction * dim(log)[1]) + 1

# Extract MCC trees
# Path to local TreeAnnotator executable
treeAnnotatorPath <- "/path/to/TreeAnnotator"

mcc_search_dir <- getwd()  # absolute path
files <- list.files(mcc_search_dir, recursive = TRUE, full.names = TRUE)

trees <- files[grepl("EmpTrees", basename(files)) & grepl("\\.trees$", basename(files))]

for (t in trees) {
  out <- sub("\\.trees$", "_MCC.tree", t)
  
  cmd <- paste0("\"", treeAnnotatorPath, "\" ","-burninTrees ", burnIn, " -heights keep ",
                "\"", t, "\" ","\"", out, "\"")
  
  message("Running: ", cmd)
  system(cmd)
}

# Extracting spatio-temporal information embedded in MCC and posterior trees of the RRW analysis
Treefile <- list.files(paste0(DirProvDTA, DTA_analysis), pattern = "EmpTrees.trees")

# Load trait and location data for two main NL Clades
df1 <- read_delim(paste0(DirProvDTA, "ProvincesDTA.NLClade1.txt"), delim = "\t")
df2 <- read_delim(paste0(DirProvDTA, "ProvincesDTA.NLClade2.txt"), delim = "\t")

NLClades <- list(
  Clade1 = df1,
  Clade2 = df2
)

for (i in 1:length(Treefile)){
  Tree_name <- sub("^(Clade[0-9]+)_ProvincesDTA.*\\.trees$", "\\1", Treefile[i])
  localTreesDirectory <- paste0(DirProvDTA,DTA_analysis,Tree_name, "_ext")
  dir.create(localTreesDirectory,showWarnings=F)
  trees = scan(paste0(DirProvDTA,DTA_analysis,Treefile[i]), what="", sep="\n", quiet=T, blank.lines.skip=F)
  index1 = which(trees=="\t\t;")[length(which(trees=="\t\t;"))]
  index2 = index1 + burnIn + 1
  indices3 = which(grepl("tree STATE", trees)); index3 = indices3[length(indices3)]
  interval = floor((index3-(index1+burnIn))/nberOfTreesToSample)
  indices = seq(index3-((nberOfTreesToSample-1)*interval),index3,interval)
  selected_trees = c(trees[c(1:index1,indices)],"End;")
  write(selected_trees, paste0(DirProvDTA,DTA_analysis,Tree_name, "_ext/", Tree_name, "_selected900.trees"))
  
  # load the subsampled trees with annotations
  subsampled_tree <- readAnnotatedNexus(paste0(DirProvDTA,DTA_analysis,Tree_name,
                                               "_ext/", Tree_name, "_selected900.trees"))
  
  for (j in 1:length(subsampled_tree)){
    tree = subsampled_tree[[j]]
    mostRecentSamplingDatum <- max(
      lubridate::decimal_date(as.Date(unlist(NLClades[[Tree_name]][, "date"]))),
      na.rm = TRUE
    )
    tab = DTA_tree_extraction1(tree, mostRecentSamplingDatum)
    tab$cladeID = rep(Tree_name, dim(tab)[1])
    write.csv(tab, paste0(localTreesDirectory,"/TreeExtractions_",j,".csv"), row.names=F, quote=F)
  }
}

# create transition matrix
Provinces <- unique(unlist(lapply(NLClades, function(df) df$location)))
Provinces

for (i in 1:length(Treefile)){
  Tree_name <- sub("^(Clade[0-9]+)_ProvincesDTA.*\\.trees$", "\\1", Treefile[i])
  localTreesDirectory <- paste0(DirProvDTA,DTA_analysis,Tree_name, "_ext")
  matrices = list()
  for (j in 1:nberOfExtractionFiles){
    mat = matrix(0, nrow=length(Provinces), ncol=length(Provinces))
    row.names(mat) = Provinces; colnames(mat) = Provinces
    tab = read.csv(paste0(localTreesDirectory, "/TreeExtractions_",j,".csv"), head=T)
    for (k in 1:nrow(tab)){
      index1 = which(Provinces==tab[k,"startLoc"])
      index2 = which(Provinces==tab[k,"endLoc"])
      mat[index1,index2] = mat[index1,index2]+1
    }
    matrices[[j]] = mat
  }
  saveRDS(matrices, paste0(DirProvDTA,DTA_analysis,Tree_name,"_matrices.rds"))
}

data <- readRDS(paste0(DirProvDTA,DTA_analysis,Tree_name,"_matrices.rds"))  # Load the RDS file (just to check)
print(data)  # View the contents

# create Bayes Factor table for each transition for the DTA and TSW analysis
for (i in 1:length(Treefile)){
  Tree_name <- sub("^(Clade[0-9]+)_ProvincesDTA.*\\.trees$", "\\1", Treefile[i])
  
  log_DTA = scan(paste0(DirProvDTA,DTA_analysis,Tree_name,"_ProvincesDTA_EmpTrees.log"), what="", sep="\n", quiet=T, blank.lines.skip=F)
  index1 = 4+burnIn; index2 = length(log_DTA); interval = round((index2-index1)/nberOfTreesToSample)
  indices = seq(index2-((nberOfTreesToSample-1)*interval),index2,interval)
  write(log_DTA[c(4,indices)],file = paste0(DirProvDTA,DTA_analysis,Tree_name,"_selected900.log"))
  log_DTA <- read_tsv(paste0(DirProvDTA,DTA_analysis,Tree_name,"_selected900.log"))
  
  log_TSW <- scan(paste0(DirProvDTA,TSW_analysis,Tree_name,"_ProvincesDTA_EmpTrees_TSW.log"), what="", sep="\n", quiet=T, blank.lines.skip=F)
  index1 = 4+burnIn; index2 = length(log_TSW); interval = round((index2-index1)/nberOfTreesToSample)
  indices = seq(index2-((nberOfTreesToSample-1)*interval),index2,interval)
  write(log_TSW[c(4,indices)],file = paste0(DirProvDTA,TSW_analysis,Tree_name,"_selected900.log"))
  log_TSW <- read_tsv(paste0(DirProvDTA,TSW_analysis,Tree_name,"_selected900.log"))
  
  BF_DTA = matrix(nrow=length(Provinces), ncol=length(Provinces))
  BF_TSW = matrix(nrow=length(Provinces), ncol=length(Provinces))
  row.names(BF_DTA) = Provinces; colnames(BF_DTA) = Provinces
  row.names(BF_TSW) = Provinces; colnames(BF_TSW) = Provinces
  for (i in 1:length(Provinces)){
    for (j in 1:length(Provinces)){
      if (i != j){
        colName = paste0("location.indicators.",gsub(" ",".",Provinces[i]),".",gsub(" ",".",Provinces[j]))
        index1 = which(colnames(log_DTA)==colName); index2 = which(colnames(log_TSW)==colName)
        p = sum(log_DTA[,index1]==1)/dim(log_DTA)[1]
        K = length(Provinces)
        q = (log(2)+K-1)/(K*(K-1))
        BF_DTA[i,j] = (p/(1-p))/(q/(1-q))
        p1 = sum(log_DTA[,index1]==1)/dim(log_DTA)[1]
        p2 = sum(log_TSW[,index2]==1)/dim(log_TSW)[1]
        BF_TSW[i,j] = (p1/(1-p1))/(p2/(1-p2))
      }
    }
  }
  write.table(round(BF_DTA,1), paste0(DirProvDTA,DTA_analysis,Tree_name,"_BF_values.csv"), sep=",", quote=F)
  write.table(round(BF_TSW,1), paste0(DirProvDTA,TSW_analysis,Tree_name,"_BF_values.csv"), sep=",", quote=F)
}

#extract info from MCC trees from DTA analysis
for (i in 1:length(Treefile)){
  Tree_name <- sub("^(Clade[0-9]+)_ProvincesDTA.*\\.trees$", "\\1", Treefile[i])
  mostRecentSamplingDatum <- max(
    lubridate::decimal_date(as.Date(unlist(NLClades[[Tree_name]][, "date"]))),
    na.rm = TRUE
  )  
  mcc_tre <- readAnnotatedNexus(paste0(DirProvDTA,DTA_analysis,Tree_name,"_ProvincesDTA_EmpTrees_MCC.tree"))
  mcc_tab<- DTA_tree_extraction1(mcc_tre, mostRecentSamplingDatum)
  mcc_tab$cladeID = rep(Tree_name, dim(mcc_tab)[1])
  write.csv(mcc_tab, paste0(DirProvDTA,DTA_analysis,Tree_name,".csv"), row.names=F, quote=F)
}

# Four-panel figure: standard DTA and adjusted-BF results for both clades
admin2 <- shapefile("Shapefiles/provincie_gegeneraliseerd.shp")
centroids = raster::coordinates(admin2) #select the centroid position of each polygon
row.names(centroids) = gsub(" ","",admin2@data$statcode)
centroids = centroids[Provinces,] #order!

onlyInternalNodesOfTipBranches = FALSE

CladeNames <- c("Clade1", "Clade2")
matrices_list <- list()
all_diag_vals <- c()
all_vals <- c()

# Set up PDF for 4-panel output
pdf(paste0(DirProvDTA, "Figure_DTA_4panels.pdf"), width=8, height=8)
par(mfrow=c(2,2), mar=c(0,0,0,0),oma=c(0,0,0,0), xaxs = "i", yaxs = "i", bty = "n", col="gray30")  # 2 rows, 2 columns layout

# Load and compute global min/max for scaling
mat1 <- read_rds(paste0(DirProvDTA,DTA_analysis,"Clade1","_matrices.rds")) 
mat2 <- read_rds(paste0(DirProvDTA,DTA_analysis,"Clade2","_matrices.rds"))

all_matrices <- list(mat1, mat2)

mean_mats <- lapply(all_matrices, function(mats) Reduce(`+`, mats) / length(mats))

# Compute global scaling parameters
global_min_diag <- min(sapply(mean_mats, function(m) min(diag(m))))
global_max_diag <- max(sapply(mean_mats, function(m) max(diag(m))))
global_min <- min(sapply(mean_mats, function(m) min(m, na.rm = TRUE)))
global_max <- max(sapply(mean_mats, function(m) max(m, na.rm = TRUE)))

legend_arrow_values <- c(1, 5, 10)
legend_cex_values   <- c(10, 50, 100) 

# Plotting parameters
multiplier1 <- 60  # point size
multiplier2 <- 20  # arrow width
multiplier3 <- 0.6 # arrow length

adjustedBFs_option <- c(FALSE, TRUE)
for (i in seq_along(CladeNames)) {
  for (adjustedBF in adjustedBFs_option) {
    
    clade <- CladeNames[i]
    matrices <- all_matrices[[i]]
    matrix_mean <- round(Reduce(`+`, matrices) / length(matrices), 1)
    mat <- matrix_mean
    
    # Read corresponding BF file
    BF_file <- if (adjustedBF) {
      paste0(DirProvDTA, TSW_analysis, clade, "_BF_values.csv")
    } else {
      paste0(DirProvDTA, DTA_analysis, clade, "_BF_values.csv")
    }
    BFs <- read.csv(BF_file, header=TRUE)
    
    # Base map
    eval(BaseMapLab)
    
    # Points for self-transitions (diagonal)
    points(
      centroids,
      cex = sqrt((multiplier1 * ((diag(mat) - global_min_diag) / (global_max_diag - global_min_diag))) / pi),
      pch = 16, col = "#cc1337"
    )
    
    # Arrows for transitions
    for (j in 1:nrow(admin2)) {
      for (k in 1:nrow(admin2)) {
        if ((j != k) && (mat[j, k] >= 1) && (!is.na(BFs[j, k])) && (BFs[j, k] > 3)) {
          LWD <- (((mat[j, k] - global_min) / (global_max - global_min)) * multiplier2) + 0.1
          arrow_len <- (multiplier3 * (mat[j, k] / global_max)) + 0.04
          curvedarrow(
            centroids[j,], centroids[k,],
            arr.length = arrow_len * 1.3, arr.width = arrow_len,
            lwd = LWD, lty = 1, lcol = "gray40", arr.col = "gray40",
            arr.pos = 0.5, curve = 0.25, dr = NA, endhead = FALSE,
            arr.type = "triangle"
          )
        }
      }
    }
    
    # Title
    mtext(paste(clade, ifelse(adjustedBF, ": DTA with adjusted BF", ": DTA")),
          side=3, line=0, cex=0.7, font=2)
  }
  # Legend
  if (clade == "Clade2" && adjustedBF) {
    
    points(cbind(rep(4.4, 4), rep(51.12, 4)),
           cex = sqrt((multiplier1 * ((legend_cex_values - global_min_diag) / (global_max_diag - global_min_diag))) / pi),
           pch = 16, col = "#cc1337", lwd = 0.3)
    
    text(4.4, 51.05, legend_cex_values[1], cex = 0.6, pos = 4)
    text(4.4, 51.1, legend_cex_values[2], cex = 0.6, pos = 4)
    text(4.4, 51.15, legend_cex_values[3], cex = 0.6, pos = 4)
    text(4.4, 51.2, legend_cex_values[4], cex = 0.6, pos = 4)
    
    # Arrows (between-location transitions)
    for (k in seq_along(legend_arrow_values)) {
      vS <- legend_arrow_values[k]
      LWD <- (((vS - global_min) / (global_max - global_min)) * multiplier2) + 0.1
      arrow <- (multiplier3 * (vS / global_max)) + 0.04
      curvedarrow(cbind(4.20, 51.0 - 0.05 * (k - 1)),
                  cbind(4.60, 51.0 - 0.05 * (k - 1)),
                  arr.length = arrow * 1.3, arr.width = arrow, lwd = LWD, lty = 1,
                  lcol = "gray40", arr.col = "gray40", arr.pos = 0.52, curve = 0,
                  dr = NA, endhead = FALSE, arr.type = "triangle")
      text(4.2, 51.0 - 0.05 * (k - 1), vS, cex = 0.6, pos = 4)
    }
  } 
}
dev.off()

## Fig. 5: Environmental rasters --------------------------------------------

# Area within which we extracted the environmental data 
NLadm1 <- st_read("Shapefiles/provincie_gegeneraliseerd.shp")
# Get bounding box
bbox <- st_bbox(NLadm1 )
# Expand the bbox 
exp_bbox <- bbox
exp_bbox["xmin"] <- bbox["xmin"] - 2.25
exp_bbox["xmax"] <- bbox["xmax"] + 2.25
exp_bbox["ymin"] <- bbox["ymin"] - 0.5
exp_bbox["ymax"] <- bbox["ymax"] + 0.5
exp_bbox
extent <- st_as_sfc(st_bbox(exp_bbox), crs = 4326)

# Prepare outline European boundaries
EUadm0 <- st_read("Shapefiles/CNTR_RG_01M_2020_4326_NLAdapt.shp")
EUadm0.1 <- crop(EUadm0, exp_bbox)
EUadm0.1 <- mask(EUadm0.1, exp_bbox)
crs(EUadm0)

plot(EUadm0.1, col="#F2F2F2", border=NA, lwd=0.1, add=F)
dev.off()

EnvRasterDir <- file.path(base_dir, "EnvRasters", "StudyRasters")

envVariableFiles = c("Urban","Agricultural","OpenVegetation","Forests",
                     "InlandWetlands","InlandWaters","GDD","MinWinterT","LogPopulationDensity")

envVariableNames = c("Urban","Agricultural","Open vegetation","Forests",
                      "Inland wetlands","Inland waters","N growing degree-days","Min winter temperature","Human pop density (log)")

rS = list(); cols = list(); colour1 = "gray98"
for (i in 1:length(envVariableFiles)){
  rS[[i]] = raster(paste0(EnvRasterDir,envVariableFiles[i],".tif"))
}

# Define colour scales
cols[[1]] = c("#FAFAFA", colorRampPalette(brewer.pal(9, "PuRd"))(99))      # Urban
cols[[2]] <- c("#FAFAFA",colorRampPalette(brewer.pal(9, "YlOrBr")[1:7])(100))    # Agricultural
cols[[3]] = c("#FAFAFA", colorRampPalette(brewer.pal(9, "Greens"))(100))    # OpenVegetation
cols[[4]] = c("#FAFAFA", colorRampPalette(brewer.pal(9, "BuGn"))(100))       # Forests
cols[[5]] = c("#FAFAFA", colorRampPalette(brewer.pal(9, "GnBu"))(100))       # InlandWetlands
cols[[6]] = c("#FAFAFA", colorRampPalette(brewer.pal(9, "Blues"))(100))    # InlandWaters
cols[[7]] = colorRampPalette(brewer.pal(9, "YlOrRd"))(100)   # GDD
cols[[8]] <- rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(100))   # MinWinterTemp (cold)
cols[[9]] = c("#FAFAFA", colorRampPalette(brewer.pal(9, "BuPu"))(99))  # PopDensity (log)

# 3×3 panel layout, one panel per environmental variable;
pdf(paste0(FigRepo,"Fig5_EnvFactors.pdf"),width=8.27, height=5.83); colNA = "grey90"
par(mfrow=c(3,3), mar=c(0,0,0,0), oma=c(1,1,1.5,1), mgp=c(0,0.4,0), lwd=0.2, bty="o", col="gray30")
for (i in 1:length(rS)){
  plot(extent, col=NA, border=NA, lwd=0.001)
  plot(rS[[i]], bty="n", box=F, axes=F, legend=F, col=cols[[i]], colNA=colNA, add=T)
  plot(rS[[i]], legend.only=T, add=T, col=cols[[i]], legend.width=0.5, legend.shrink=0.3, smallplot=c(0.936,0.951,0.035,0.965),
       legend.args=list(text="", cex=0.7, line=0.3, col="gray30"), horizontal=F,
       axis.args=list(cex.axis=0.75, lwd=0, lwd.tick=0.2, col.tick="gray30", tck=-1.0, col="gray30", col.axis="gray30", line=0, mgp=c(0,0.50,0)))
  plot(EUadm0.1, col=NA, border="#737373", lwd=0.4, add=T)
  plot(extent, col=NA, border="#737373", lwd=0.001, add=T)
  usr <- par("usr")
  mtext(envVariableNames[i],
        side = 3,
        at   = mean(usr[1:2]),  # x-center of current panel
        line = -1.6,            
        cex  = 0.7,
        font = 1, 
        col  = "#737373")}
dev.off()
