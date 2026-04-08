####################################
#   Mutation richness analysis     #
#         Version 1.0              #
# Author: Adam A. Capoferri, PhD   #
# Contact: adam.capoferri@nih.gov  #
####################################

# This script is a sliding window analysis to illustrate where mutations of choice (in this case) GR (R=G/A) positions are in the genome. This is meant to help indicate regions that might be targeted through APOBEC3 activity (GR>AR)

# load libraries
library(Biostrings)
library(GenomicRanges)

# load genome of interest, change as necessary
genome <- readDNAStringSet("filename.fasta")

# set the window and step size
scan_windows <- function(seq, window_size = 300, step = 10) { 
  
  starts <- seq(1, length(seq) - window_size + 1, by = step)
  
  res <- data.frame(
    start = starts,
    end = starts + window_size - 1,
    G_count = NA,
    GR_count = NA,
    GR_fraction = NA
  )
  
  for (i in seq_along(starts)) {
    w <- subseq(seq, starts[i], starts[i] + window_size - 1)
    
    g_count <- countPattern("G", w)
    
    gr_count <- countPattern("GA", w) + countPattern("GG", w)
    
    res$G_count[i] <- g_count
    res$GR_count[i] <- gr_count
    res$GR_fraction[i] <- gr_count / length(w)
  }
  
  res
}

results <- lapply(genome, scan_windows)
results <- do.call(rbind, results)

find_GR_sites <- function(seq) {
  
  ga <- matchPattern("GA", seq)
  gg <- matchPattern("GG", seq)
  
  c(start(ga), start(gg))
}

gr_sites <- lapply(genome, find_GR_sites)

plot(results$start, results$GR_fraction, type="l",
     xlab="Genome position",
     ylab="GR fraction",
     main="GR richness across genome"),
    yaxis

high_GR <- results$GR_fraction > 0.10
points(results$start[high_GR], results$GR_fraction[high_GR], col="red", pch=16)

#### END ####
