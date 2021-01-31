# FROM CLUSTAL TO HEATMAP
# This script uses as input an alignment file in Clustal format,
# as the one produced by default by EBI Clustal Omega:
# https://www.ebi.ac.uk/Tools/msa/clustalo/
# The alignment is splitted in genomic windows with size chosen
# by the user, the percentage of identical characters between a
# reference sequence is calculated, and a heatmap is plotted.
# Moreno Colaiacovo - 2021

# PACKAGES
library(seqinr) # To read the alignment file
library(ggplot2) # To draw the heatmap
library(showtext) # To use additional fonts in the plot
library(viridis) # To use the color palette

# PLOT FUNCTION
# Parameters:
# AL: alignment object
# FILE_NAME: name of the output file
# WS: size of the genomic window
# nS: nucleotide start for the heatamp (useful to zoom in specific regions)
# nE: nucleotide end for the heatmap (useful to zoom in specific regions)
# S: index of the reference sequence in AL$nam
# seq_order: list with sequence names taken from AL$nam in the preferred order
# min_color: identity % corresponding to the lower value in the color scale
# plot_break: nucleotide distance between ticks in x-axis
# plot_title: title of the plot
plot_heatmap <- function(AL, FILE_NAME, WS, nS, nE, S, seq_order, min_color, plot_break, plot_title) {
  
  df <- data.frame(POS=integer(), IDENTITY=double(), NAME=character())

  for (V in 1:length(AL$nam)) {
    for (n in seq(nS,nE,WS)) {
      S1 <- substr(AL$seq[S],n,n+WS-1)
      S2 <- substr(AL$seq[V],n,n+WS-1)
      ids <- 0
      for (x in 1:nchar(S1)) {
        if (substr(S1,x,x) == substr(S2,x,x)) {ids <- ids + 1}
      }
      ID <- 100*(ids/WS)
      df <- rbind(df,data.frame(POS=n,IDENTITY=ID,NAME=AL$nam[V]))
    }
  }

  P <- ggplot(data=df[df$NAME!=AL$nam[S],], aes(POS, factor(NAME,levels=rev(reordered)), fill= IDENTITY)) +
    geom_tile(colour = "white", size=2) +
    ggtitle(plot_title) +
    scale_fill_viridis(option="D",na.value='white',limits=c(min_color,100)) +
    labs(x="",y="") +
    scale_x_continuous(position='top', breaks = seq(0, 30000, by = plot_break)) + 
    theme_classic() +
    theme(panel.grid = element_blank(),
          plot.margin=unit(c(1,1,1,1),"cm"),
          legend.title=element_blank(),
          legend.position="right",
          legend.margin=margin(grid::unit(0,"cm")),
          legend.text=element_text(size=100),
          legend.key.height=grid::unit(6,"cm"),
          legend.key.width=grid::unit(2,"cm"),
          axis.title=element_text(size=100,family="Roboto",margin = unit(c(4, 0, 4, 0), "cm")),
          plot.title = element_text(size=150,family='Roboto'),
          axis.text.x = element_text(size=100,family='Roboto',angle=45, hjust=0),
          axis.text.y = element_text(size=100,family='Roboto', margin = unit(c(0, 0.5, 0, 0), "cm")),
          axis.ticks.length = unit(0.5, "cm"))

  ggsave(FILE_NAME,P,dev="png", width = 35, height = 20, units='in', bg='transparent')

}

# MAIN PROGRAM

# ACTIVATE SHOWTEXT
showtext_auto()
x11()

# LOAD FONT
# The Roboto-Light.ttf can be downloaded from Google Fonts:
# https://fonts.google.com/specimen/Roboto
# Remember to use the correct path when calling font_add
font_add(family='Roboto', regular='Roboto-Light.ttf')

# SET WORKING DIRECTORY
# Set your working directory here
setwd("~/GIT_folder/from_Clustal_to_Heatmap")

# LOAD ALIGNMENT
# The Clustal alignment should be placed in the working directory
AL <- read.alignment('SARS2_relatives.aln', 'clustal')

# REMOVE TABS
AL$seq <- lapply(AL$seq, function(y) gsub("\t", "", y,fixed=TRUE))

# GET INDEX OF REFERENCE SEQUENCE
# Choose a reference sequence in the alignment and extract its index in AL$nam
S <- match('SARS2',AL$nam)

# ORDER VIRUSES
# Choose the order of the sequences in the heatmap.
reordered = c('RaTG13','RShSTT182','RShSTT200',
              'RacCS203','PnMP789','ZC45','ZXC21',
              'PnGX-P2V_2018','PnGX-P1E_2017')


# PLOT HEATMAP FOR DIFFERENT GENOMIC REGIONS
plot_heatmap(AL,'SARS2_relatives.Whole_Genome.png',750, 1, nchar(AL$seq[S]), S, reordered, 50, 2000, 'Genome identity of SARS-CoV-2 relatives - Whole genome')
plot_heatmap(AL,'SARS2_relatives.S_gene.png',100, 21635, 25490, S, reordered, 0, 500, 'Genome identity of SARS-CoV-2 relatives - S gene')
plot_heatmap(AL,'SARS2_relatives.RBD.png',25, 22575, 23417, S, reordered, 0, 100, 'Genome identity of SARS-CoV-2 relatives - Receptor Binding Domain')
plot_heatmap(AL,'SARS2_relatives.RBM.png',5, 22977, 23192, S, reordered, 0, 10, 'Genome identity of SARS-CoV-2 relatives - Receptor Binding Motif')

# TURN OFF SHOWTEXT
showtext_auto(FALSE)

### GENOMIC COORDINATES ###
# The following are the genomic coordinates for
# S gene, RBD and RBM. Please note that the coordinates
# to be used in the plot_heatmap function are different,
# because they refer to the alignment positions.
# S gene
#21563-25384 in SARS2 genome
# RBD
#22469-23311 in SARS2 genome
# RBM
#22871-23086 in SARS2 genome
