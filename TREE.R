#if (!requireNamespace("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")

#BiocManager::install("Biostrings")


#install.packages("ape")
#install.packages("ggplot2")
#install.packages("ggtree")
#install.packages("phangorn")
#install.packages("seqinr")

library("ape")
library("Biostrings")
library("ggplot2")
library("ggtree")
library("phangorn")
library("seqinr")

#build ![tree](https://a-little-book-of-r-for-bioinformatics.readthedocs.io/en/latest/src/chapter5.html)

coral_co1 <- read.alignment("CO1_DNA.phy", "phylip")

#creates the function to clean the alignment data
cleanAlignment <- function(alignment, minpcnongap, minpcid)
{
  # make a copy of the alignment to store the new alignment in:
  newalignment <- alignment
  # find the number of sequences in the alignment
  numseqs <- alignment$nb
  # empty the alignment in "newalignment")
  for (j in 1:numseqs) { newalignment$seq[[j]] <- "" }
  # find the length of the alignment
  alignmentlen <- nchar(alignment$seq[[1]])
  # look at each column of the alignment in turn:
  for (i in 1:alignmentlen)
  {
    # see what percent of the letters in this column are non-gaps:
    nongap <- 0
    for (j in 1:numseqs)
    {
      seqj <- alignment$seq[[j]]
      letterij <- substr(seqj,i,i)
      if (letterij != "-") { nongap <- nongap + 1}
    }
    pcnongap <- (nongap*100)/numseqs
    # Only consider this column if at least minpcnongap % of the letters are not gaps:
    if (pcnongap >= minpcnongap)
    {
      # see what percent of the pairs of letters in this column are identical:
      numpairs <- 0; numid <- 0
      # find the letters in all of the sequences in this column:
      for (j in 1:(numseqs-1))
      {
        seqj <- alignment$seq[[j]]
        letterij <- substr(seqj,i,i)
        for (k in (j+1):numseqs)
        {
          seqk <- alignment$seq[[k]]
          letterkj <- substr(seqk,i,i)
          if (letterij != "-" && letterkj != "-")
          {
            numpairs <- numpairs + 1
            if (letterij == letterkj) { numid <- numid + 1}
          }
        }
      }
      pcid <- (numid*100)/(numpairs)
      # Only consider this column if at least %minpcid of the pairs of letters are identical:
      if (pcid >= minpcid)
      {
        for (j in 1:numseqs)
        {
          seqj <- alignment$seq[[j]]
          letterij <- substr(seqj,i,i)
          newalignmentj <- newalignment$seq[[j]]
          newalignmentj <- paste(newalignmentj,letterij,sep="")
          newalignment$seq[[j]] <- newalignmentj
        }
      }
    }
  }
  return(newalignment)
}

cleanedcoral_co1 <- cleanAlignment(coral_co1,30,30)

printMultipleAlignment <- function(alignment, chunksize=60)
{
  # this function requires the Biostrings package
  require("Biostrings")
  # find the number of sequences in the alignment
  numseqs <- alignment$nb
  # find the length of the alignment
  alignmentlen <- nchar(alignment$seq[[1]])
  starts <- seq(1, alignmentlen, by=chunksize)
  n <- length(starts)
  # get the alignment for each of the sequences:
  aln <- vector()
  lettersprinted <- vector()
  for (j in 1:numseqs)
  {
    alignmentj <- alignment$seq[[j]]
    aln[j] <- alignmentj
    lettersprinted[j] <- 0
  }
  # print out the alignment in blocks of 'chunksize' columns:
  for (i in 1:n) { # for each of n chunks
    for (j in 1:numseqs)
    {
      alnj <- aln[j]
      chunkseqjaln <- substring(alnj, starts[i], starts[i]+chunksize-1)
      chunkseqjaln <- toupper(chunkseqjaln)
      # Find out how many gaps there are in chunkseqjaln:
      gapsj <- countPattern("-",chunkseqjaln) # countPattern() is from Biostrings package
      # Calculate how many residues of the first sequence we have printed so far in the alignment:
      lettersprinted[j] <- lettersprinted[j] + chunksize - gapsj
      print(paste(chunkseqjaln,lettersprinted[j]))
    }
    print(paste(' '))
  }
}

printMultipleAlignment(cleanedcoral_co1)

coral_co1bin <- as.DNAbin(coral_co1)
coral_co1Dist <- dist.dna(coral_co1bin)
coral_co1Dist

rootedNJtree <- function(alignment, theoutgroup, type)
{
  # load the ape and seqinR packages:
  require("ape")
  require("seqinr")
  # define a function for making a tree:
  makemytree <- function(alignmentmat, outgroup=`theoutgroup`)
  {
    alignment <- ape::as.alignment(alignmentmat)
    if      (type == "protein")
    {
      mydist <- dist.alignment(alignment)
    }
    else if (type == "DNA")
    {
      alignmentbin <- as.DNAbin(alignment)
      mydist <- dist.dna(alignmentbin)
    }
    mytree <- nj(mydist)
    mytree <- makeLabel(mytree, space="") # get rid of spaces in tip names.
    myrootedtree <- root(mytree, outgroup, r=TRUE)
    return(myrootedtree)
  }
  # infer a tree
  mymat  <- as.matrix.alignment(alignment)
  myrootedtree <- makemytree(mymat, outgroup=theoutgroup)
  # bootstrap the tree
  myboot <- boot.phylo(myrootedtree, mymat, makemytree)
  # plot the tree:
  plot.phylo(myrootedtree, type="p")  # plot the rooted phylogenetic tree
  nodelabels(myboot,cex=0.7)          # plot the bootstrap values
  mytree$node.label <- myboot   # make the bootstrap values be the node labels
  return(mytree)
}

unrootedNJtree <- function(alignment,type)
{
  # this function requires the ape and seqinR packages:
  require("ape")
  require("seqinr")
  # define a function for making a tree:
  makemytree <- function(alignmentmat)
  {
    alignment <- ape::as.alignment(alignmentmat)
    if      (type == "protein")
    {
      mydist <- dist.alignment(alignment)
    }
    else if (type == "DNA")
    {
      alignmentbin <- as.DNAbin(alignment)
      mydist <- dist.dna(alignmentbin)
    }
    mytree <- nj(mydist)
    mytree <- makeLabel(mytree, space="") # get rid of spaces in tip names.
    return(mytree)
  }
  # infer a tree
  mymat  <- as.matrix.alignment(alignment)
  mytree <- makemytree(mymat)
  # bootstrap the tree
  myboot <- boot.phylo(mytree, mymat, makemytree)
  # plot the tree:
  plot.phylo(mytree,type="u")   # plot the unrooted phylogenetic tree
  nodelabels(myboot,cex=0.7)    # plot the bootstrap values
  mytree$node.label <- myboot   # make the bootstrap values be the node labels
  return(mytree)
}

coral_co1Tree <- unrootedNJtree(cleanedcoral_co1, type = "DNA")
write.tree(coral_co1Tree, "Cleanedcoral_co1Tree.tre")

#build ![Tree](https://www.molecularecologist.com/2016/02/quick-and-dirty-tree-building-in-r/)
#Haven't got the below to work

CO1 <- read.dna("CO1_DNA.phy", format="interleaved")
CO1_phyDat <- phyDat(CO1, type = "DNA", levels = NULL)
CO110 <- subset(CO1_phyDat, 1:10)
CO110_phyDat <- phyDat(CO110, type = "DNA", levels = NULL)

mt <- modelTest(CO110)
print(mt)
#dna_dist <- dist.ml(CO110, model="JC69")

CO1_UPGMA <- upgma(coral_co1Dist)
CO1_NJ  <- NJ(coral_co1Dist)
plot(CO1_UPGMA, main="UPGMA")
plot(CO1_NJ, main = "Neighbor Joining")

#parsimony(CO1_UPGMA)
#parsimony(CO1_NJ)
#CO1_optim <- optim.parsimony(CO1_NJ, CO1)
#CO1_pratchet <- pratchet(CO110)
#plot(CO1_optim)
#plot(CO1_pratchet)

fit <- pml(CO1_UPGMA, CO1_phyDat)
print(fit)
fitJC <- optim.pml(fit, model = "JC", rearrangement = "stochastic")
logLik(fitJC)
bs <- bootstrap.pml(fitJC, bs=100, optNni=TRUE, multicore=TRUE, control = pml.control(trace=0))
plotBS(midpoint(fitJC$tree), bs, p = 50, type="p")
