
#wg-mining-validation-source.R
#
# This script provides the source code necessary to analyses the metabarcoding and hts data from the faecal, water and sediment samples collected in Slimbridge between 2017-2020. 
# Copyright (c) Sarah Nichols, 2022, except where indicated
# Date Created: 2022-09-7


# --------------------------------------------------------------------------
# REQUIRES
# --------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(devtools)
  library(ggplot2)
  library(dada2)
  library(dplyr)
  library(purrr)
  library(phyloseq)
  library(patchwork)
  library(vegan)
  library(eulerr)
  library(sf)
  library(sp)
  library(ggspatial)
  library(Biostrings); packageVersion("Biostrings")
  library(ShortRead); packageVersion("ShortRead")
  library(reshape2); packageVersion("reshape2")
  library(gridExtra); packageVersion("gridExtra")
  library(paletteer)
  library(scales)
  library(colorspace)
  library(ggrepel)
  library(caret)
  library(readr)
  library(ggpubr)
  library(taxizedb)
  library(ComplexHeatmap)
  library(grid)

  })



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

text_size=15

theme_set(theme_bw())

scale_colour_discrete <-  function(palname=pal, ...){
  scale_colour_brewer(palette=palname, ...)
}
scale_fill_discrete <-  function(palname=pal, ...){
  scale_fill_brewer(palette=palname, ...)
}

palette <- paletteer_d(palette = "khroma::stratigraphy", n = 175)


methods_pal <- c(palette[c(168,74,64, 133, 18)])
#host_pal <- c(palette[c(102,138,109, 146)])
#trans_pal <- c(palette[c(3, 21, 82, 87)])

palette_2 <- paletteer::paletteer_d("ggthemes::Tableau_20")

host_pal <- palette_2[c(5, 6, 7, 8)]
trans_pal <- palette_2[c(1, 2, 11, 12)]


# Function to generate shades of a base color, avoiding white
generate_shades_lighten <- function(base_color, num_shades) {
  # Validate the input color
  if (!grepl("^#", base_color)) {
    stop("Invalid base color. Please provide a hex color code (e.g., '#RRGGBB').")
  }
  
  # Generate a sequence of lightness values, avoiding extremes (0.9 prevents white)
  lightness_steps <- seq(0.2, 0.7, length.out = num_shades)  # Adjust range as needed
  
  # Apply lightening to the base color
  shades <- sapply(lightness_steps, function(lightness) {
    lighten(base_color, amount = lightness)
  })
  
  # Return the generated shades
  return(shades)
}



# Function to generate shades of a base color, avoiding white
generate_shades_darken <- function(base_color, num_shades) {
  # Validate the input color
  if (!grepl("^#", base_color)) {
    stop("Invalid base color. Please provide a hex color code (e.g., '#RRGGBB').")
  }
  
  # Generate a sequence of lightness values, avoiding extremes (0.9 prevents white)
  darkness_steps <- seq(0.2, 0.7, length.out = num_shades)  # Adjust range as needed
  
  # Apply lightening to the base color
  shades <- sapply(darkness_steps, function(lightness) {
    darken(base_color, amount = lightness)
  })
  
  # Return the generated shades
  return(shades)
}


# Generate 5 shades of the base color
#unassembled_pal <- generate_shades_lighten(c(palette[74]), 3)
#assembled_pal <- generate_shades_lighten(c(palette[168]), 3)

unassembled_pal <- generate_shades_lighten("#79706EFF", 3)
assembled_pal <- generate_shades_lighten(c("#79706EFF"), 3)

#heatmap_pal <- c(palette[c(44,20,30, 85)])
heatmap_pal <- (c(rev(generate_shades_lighten(c(palette_2[17]), 2)), (generate_shades_darken(c(palette_2[17]), 2))))
total_pal <- c(heatmap_pal, methods_pal, host_pal, trans_pal)


parasite_genera_present <- c("Leishmania", "Plasmodium", "Coccidioides", "Toxoplasma", "Cryptosporidium",
                             "Babesia", "Eimeria", "Ascogregarina", "Crithidia", "Porospora", "Trypanosoma", 
                             "Cystoisospora", "Fusarium", "Phytophthora", "Theileria", "Gregarina", "Spironucleus", 
                             "Besnoitia", "Hammondia", "Cytauxzoon", "Trichomonas", "Strigomonas", "Encephalitozoon", "Cyclospora",
                             "Neospora", "Hepatocystis","Acanthamoeba", "Aspergillus", "Sarcocystis", "Amauroascus", "Entamoeba", 
                             "Hamiltosporidium", "Leptomonas", "Giardia", "Cardiosporidium", "Haemoproteus", "Nosema", 
                             "Tritrichomonas", "Enterocytozoon", "Fonsecaea", "Exophiala", "Porcisia", "Lotmaria", "Pyricularia",
                             "Endotrypanum", "Ecytonucleospora", "Sporisorium", "Vickermania", "Trypanozoon", "Hexamita",
                             "Kipferlia", "Sinuolinea", "Goussia", "Zelonia", "Parahaemoproteus", "Raphidascaris", "Proteocephalus", "Dibothriocephalus",
                             "Cystodiscus", "Rhabdias", "Philometridae", "Selenidium", "Clonorchis", "Rhytidocystis", "Platyproteum", "Isospora", "Leucocytozoon",
                             "Haemocystidium", "Blastocystis", "Nycteria", "Echinococcus", "Xiphinema", "Paratylenchus", "Olssonium", "Longidorus", "Astomonema",
                             "Tylencholaimellus", "Tetracapsuloides", "Echinobothrium", "Trichodorus", "Schistosoma", "Austramphilina", "Latyspora", "Polymorphus", 
                             "Heterodera", "Cephalenchus", "Philometra", "Ascaris", "Urogonimus", "Euclinostomum", "Meara", "Wallacemonas", "Pratylenchus", "Herpetomonas", 
                             "Novymonas", "Herpetomonas", "Jaenimonas", "Polystoma",   "Albugo" , "Hayloperonospora")



transmission_modes_dict = c(
  "Leishmania" = "flying_insect",
  "Plasmodium" = "flying_insect",
  "Coccidioides" = "direct",
  "Toxoplasma" = "two_host", 
  "Cryptosporidium" = "direct",
  "Babesia" = "tick_borne",
  "Eimeria" = "direct",
  "Ascogregarina" = "direct", 
  "Crithidia" = "direct",
  "Porospora" = "two_host",
  "Trypanosoma" = "flying_insect",
  "Cystoisospora" = "direct",
  "Fusarium" = "direct",
  "Phytophthora" = "direct",
  "Theileria" = "tick_borne",
  "Gregarina" = "direct",
  "Spironucleus" = "direct",
  "Besnoitia" = "flying_insect",
  "Hammondia" = "two_host",
  "Cytauxzoon" = "tick_borne",
  "Trichomonas" = "direct",           
  "Strigomonas" = "direct",
  "Encephalitozoon" = "two_host",                           
  "Cyclospora" = "direct",
  "Neospora" = "two_host",
  "Hepatocystis" = "flying_insect",
  "Acanthamoeba" = "direct",
  "Aspergillus" = "direct", 
  "Sarcocystis" = "two_host",
  "Amauroascus" = "direct",
  "Entamoeba" = "direct",
  "Hamiltosporidium" = "direct",
  "Leptomonas" = "direct", 
  "Giardia" = "direct",
  "Cardiosporidium" = "direct",
  "Haemoproteus" = "flying_insect",
  "Nosema" = "direct",
  "Tritrichomonas" = "direct",
  "Enterocytozoon" = "two_host",
  "Fonsecaea" = "direct",
  "Exophiala" = "direct",
  "Porcisia" = "flying_insect",
  "Lotmaria" = "direct",
  "Pyricularia" = "direct",
  "Endotrypanum" = "flying_insect",
  "Ecytonucleospora" = "direct",
  "Sporisorium" = "direct",
  "Vickermania" = "direct",
  "Trypanozoon" = "flying_insect",
  "Hexamita" = "direct",
  "Sinuolinea" = "two_host",
  "Goussia" = "direct",
  "Kipferlia" = "unknown",
  "Zelonia" = "flying_insect",
  "Parahaemoproteus" = "flying_insect",
  "Raphidascaris" = "two_host",
  "Proteocephalus" = "two_host",
  "Dibothriocephalus" = "two_host",
  "Cystodiscus" = "two_host",
  "Rhabdias" = "two_host",
  "Philometridae" = "unknown",
  "Selenidium" = "unknown",
  "Clonorchis" = "two_host",
  "Rhytidocystis" = "two_host", 
  "Platyproteum" = "two_host",
  "Isospora" = "direct or two_host",
  "Leucocytozoon" = "flying_insect",
  "Haemocystidium" = "flying_insect",
  "Blastocystis" = "direct",
  "Nycteria" = "flying_insect", 
  "Echinococcus" = "two_host",
  "Xiphinema" = "direct",
  "Paratylenchus" = "unknown", 
  "Olssonium" = "unknown",
  "Longidorus" = "unknown",
  "Astomonema" = "unknown",
  "Tylencholaimellus" = "unknown",
  "Tetracapsuloides" = "two_host",
  "Echinobothrium" = "unknown",
  "Trichodorus" = "unknown",
  "Schistosoma" = "two_host",
  "Austramphilina" = "two_host",
  "Latyspora" = "two_host",
  "Polymorphus" = "two_host", 
  "Heterodera" = "direct",
  "Cephalenchus" = "direct",
  "Philometra" = "two_host",
  "Ascaris" = "direct",
  "Urogonimus" = "two_host",
  "Euclinostomum" = "two_host",
  "Meara" = "direct",
  "Wallacemonas" = "direct",
  "Pratylenchus" = "direct",
  "Herpetomonas" = "direct",
  "Novymonas" = "direct",
  "Herpetomonas" = "direct",
  "Jaenimonas" = "direct",
  "Polystoma" = "direct",
  "Albugo" = "direct", 
  "Hayloperonospora" = "direct"
)



recorded_dict = c(
  "Leishmania" = "birds",
  "Plasmodium" = "birds",
  "Coccidioides" = "birds",
  "Toxoplasma" = "birds", 
  "Cryptosporidium" = "birds",
  "Babesia" = "birds",
  "Eimeria" = "birds",
  "Ascogregarina" = "invertebrates", 
  "Crithidia" = "invertebrates",
  "Porospora" = "invertebrates",
  "Trypanosoma" = "birds",
  "Cystoisospora" = "birds",
  "Fusarium" = "plants",
  "Phytophthora" = "plants",
  "Theileria" = "vertebrates",
  "Gregarina" = "invertebrates",
  "Spironucleus" = "birds",
  "Besnoitia" = "birds",
  "Hammondia" = "vertebrates",
  "Cytauxzoon" = "vertebrates",
  "Trichomonas" = "birds",           
  "Strigomonas" = "invertebrates",
  "Encephalitozoon" = "birds",                           
  "Cyclospora" = "birds",
  "Neospora" = "birds",
  "Hepatocystis" = "vertebrates",
  "Acanthamoeba" = "birds",
  "Aspergillus" = "birds", 
  "Sarcocystis" = "birds",
  "Amauroascus" = "vertebrates",
  "Entamoeba" = "birds",
  "Hamiltosporidium" = "invertebrates",
  "Leptomonas" = "invertebrates", 
  "Giardia" = "birds",
  "Cardiosporidium" = "birds",
  "Haemoproteus" = "birds",
  "Nosema" = "invertebrates",
  "Tritrichomonas" = "birds",
  "Enterocytozoon" = "birds",
  "Fonsecaea" = "vertebrates",
  "Exophiala" = "birds",
  "Porcisia" = "vertebrates",
  "Lotmaria" = "invertebrates",
  "Pyricularia" = "plants",
  "Endotrypanum" = "vertebrates",
  "Ecytonucleospora" = "birds",
  "Sporisorium" = "plants",
  "Vickermania" = "invertebrates",
  "Trypanozoon" = "birds",
  "Hexamita" = "birds",
  "Sinuolinea"= "vertebrates",
  "Goussia" = "vertebrates",
  "Kipferlia" = "unknown",
  "Zelonia" = "vertebrates",
  "Parahaemoproteus" = "birds",
  "Raphidascaris" = "vertebrates",
  "Proteocephalus" = "vertebrates",
  "Dibothriocephalus" = "birds",
  "Cystodiscus" = "vertebrates",
  "Rhabdias" = "vertebrates",
  "Philometridae" = "vertebrates",
  "Selenidium" = "invertebrates",
  "Clonorchis" = "vertebrates",
  "Rhytidocystis" = "vertebrates", 
  "Platyproteum" = "invertebrates",
  "Isospora" = "birds",
  "Leucocytozoon" = "birds",
  "Haemocystidium" = "vertebrates",
  "Blastocystis" = "birds",
  "Nycteria" = "vertebrates",
  "Echinococcus" = "vertebrates",
  "Xiphinema" = "plants",
  "Paratylenchus" = "plants",
  "Olssonium" = "vertebrates",
  "Longidorus" = "plants", 
  "Astomonema" = "vertebrates",
  "Tylencholaimellus" = "plants",
  "Tetracapsuloides" = "vertebrates",
  "Echinobothrium" = "vertebrates",
  "Trichodorus" = "plants",
  "Schistosoma" = "birds",
  "Austramphilina" = "vertebrates",
  "Latyspora" = "vertebrates",
  "Polymorphus" = "birds", 
  "Heterodera" = "plants",
  "Cephalenchus" = "plants",
  "Philometra" = "vertebrates", 
  "Ascaris" = "birds",
  "Urogonimus" = "birds",
  "Euclinostomum" = "birds",
  "Meara" = "invertebrates",
  "Wallacemonas" = "invertebrates",
  "Pratylenchus" = "plants",
  "Herpetomonas" = "invertebrates",
  "Novymonas" = "invertebrates",
  "Herpetomonas" = "invertebrates",
  "Jaenimonas" = "invertebrates",
  "Polystoma" = "vertebrates",
  "Albugo" = "plants", 
  "Hayloperonospora" = "plants"
)


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
