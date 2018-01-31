
# SimEli algorithm from "SimEli: Similarity Elimination Method
# for Sampling Distant Entries in Development of Core Collections"
# (Ramesh et al., 2014)
simeli <- function(dist, size, elim.crit, verbose = TRUE){
  
  # check input
  if(!is.matrix(dist) || !is.numeric(dist) || !isSymmetric(dist)){
    stop("Argument 'dist' should be a symmetric distance matrix.")
  }
  if(is.null(rownames(dist)) || is.null(colnames(dist)) || rownames(dist) != colnames(dist)){
    stop("Distance matrix should have symmetric row and column names.")
  }
  n <- nrow(dist)
  if(size < 2 || size >= n){
    stop("Core size should be at least two and less than the size of the dataset.")
  }
  
  # iterate elimination procedure
  core <- rownames(dist)
  diag(dist) <- Inf # trick to skip diagonal when looking for minimum distance
  message(sprintf("Initial size: %d", length(core)))
  while(length(core) > size){
    
    # 1: find two closest items in the current selection
    core.dist <- dist[core, core]
    min.dist <- min(core.dist)
    closest <- which(core.dist == min.dist, arr.ind = TRUE)[1,]
    item1 <- core[closest["row"]]
    item2 <- core[closest["col"]]
    elim.cand <- c(item1, item2)
      
    # 2: eliminate one item based on given elimination criterion
    elim <- elim.crit(elim.cand, core.dist)
    
    # update core
    core <- setdiff(core, elim)
    if(verbose){
      message(sprintf(
        "Closest pair: %s - %s (d = %f). Eliminate %s. New core size: %d.",
        item1, item2, min.dist, elim, length(core)
      ))
    }

  }
  
  message(sprintf("Final size: %d", length(core)))
  return(core)
  
}

# eliminate item with the smallest average distance to the remaining accessions (A-RA)
simeli.ara <- function(dist, size, verbose = TRUE, ...){
  simeli(dist, size, function(elim.cand, core.dist){
    remaining <- setdiff(rownames(core.dist), elim.cand)
    item1.A.RA <- mean(core.dist[elim.cand[1], remaining])
    item2.A.RA <- mean(core.dist[elim.cand[2], remaining])
    elim <- elim.cand[which.min(c(item1.A.RA, item2.A.RA))]
    return(elim)
  }, verbose)
}

# eliminate item with the lowest expected heterozygosity (HE)
simeli.he <- function(dist, size, verbose = TRUE, freqs, alleles){
  simeli(dist, size, function(elim.cand, core.dist){
    remaining <- setdiff(rownames(core.dist), elim.cand)
    item1.HE <- expected.heterozygosity(freqs[c(elim.cand[1], remaining),], alleles)
    item2.HE <- expected.heterozygosity(freqs[c(elim.cand[2], remaining),], alleles)
    elim <- elim.cand[which.min(c(item1.HE, item2.HE))]
    return(elim)
  }, verbose)
}

