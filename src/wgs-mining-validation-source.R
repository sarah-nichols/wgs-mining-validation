
#wg-mining-validation-source.R
#
# This script provides the source code necessary to analyses the metabarcoding and hts data from the faecal, water and sediment samples collected in Slimbridge between 2017-2020. 
# Copyright (c) Sarah Nichols, 2022, except where indicated
# Date Created: 2022-09-7


# --------------------------------------------------------------------------
# REQUIRES
# --------------------------------------------------------------------------

# suppressPackageStartupMessages({
#   library(tidyverse)
#   library(devtools)
#   library(ggplot2)
#   library(dada2)
#   library(DECIPHER)
#   library(dplyr)
#   library(purrr)
#   library(phyloseq)
#   library(RColorBrewer)
#   library(vegan)
#   library(sf)
#   library(sp)
#   library(ggspatial)
#   library(Biostrings); packageVersion("Biostrings")
#   library(ShortRead); packageVersion("ShortRead")
#   library(reshape2); packageVersion("reshape2")
#   library(gridExtra); packageVersion("gridExtra")
#   })



# --------------------------------------------------------------------------
# PATHS
# --------------------------------------------------------------------------

data_path <- file.path(getwd(), "data")
reports_path <- file.path(getwd(), "reports")
figures_path <- file.path(getwd(), "reports", "figures")

if (!dir.exists(data_path)) {
  dir.create(data_path, recursive = TRUE)
}

if (!dir.exists(reports_path)) {
  dir.create(reports_path, recursive = TRUE)
}

if (!dir.exists(figures_path)) {
  dir.create(figures_path, recursive = TRUE)
}


# --------------------------------------------------------------------------
# RESOURCES
# --------------------------------------------------------------------------

text_size=12

theme_set(theme_bw())
pal = "Set1"
scale_colour_discrete <-  function(palname=pal, ...){
  scale_colour_brewer(palette=palname, ...)
}
scale_fill_discrete <-  function(palname=pal, ...){
  scale_fill_brewer(palette=palname, ...)
}

major_taxa_colours <- c("#ABDDDE", "#78B7C5", "#3B9AB2", "#85D4E3", "#46ACC8", 
                                 "#23395D", "#FF7F94","#EFFF7F", "#FC0C60", "#E57FFF", 
                                 "#D632FF", "#8E00B2", "#F2300F","#397D7F","#CDC08C", 
                                 "#A27331", "#362204", "#64471E", "#382811", "#BFBFBF", 
                                 "#B40F20", "#E5A039", "#EAB676", "#FD5602", "#FE6E00", 
                                 "#FF8303", "#FEDEBE", "#FFAF42", "#FC9C0C", "#7FFF32", 
                                 "#42B200", "#006600", "#003314", "#42B200")
                                 
# --------------------------------------------------------------------------
# FUNCTIONS
# --------------------------------------------------------------------------

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m - sd(x)
  ymax <- m + sd(x)
  return(c(y = m, ymin = ymin, ymax = ymax))
}

'%!in%' <- function(x,y)!('%in%'(x,y))

getN <- function(x) sum(getUniques(x))

group_labeller <- function(variable,value){
  return(group_names[value])
}

gg_color_hue <- function(n){
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

display_venn <- function(x){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL)
  grid.draw(venn_object)
}

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
