library("ape")
library("Biostrings")
library("ggplot2")
library("ggtree")

nwk <- system.file("extdata", "Cleanedcoral_co1Tree.tre", package="ggtree")

tree <- read.tree(nwk)

tree
