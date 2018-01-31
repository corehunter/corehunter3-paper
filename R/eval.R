
# modified Rogers genetic distance
modified.rogers <- function(freqs1, freqs2, num.markers){
  diff <- freqs1 - freqs2
  diff[is.na(diff)] <- 0.0
  mr <- sqrt((sum(diff^2))/(2*num.markers))
  return(mr)
}

# expected proportion of heterozygous loci
expected.heterozygosity <- function(freqs, alleles){
  num.markers <- length(alleles)
  allele.counts <- sapply(alleles, length)
  freqs[is.na(freqs)] <- 0
  if(is.matrix(freqs)){
    freqs <- colMeans(freqs)
  }
  cum.allele.counts <- c(0, cumsum(allele.counts))
  f.square.sums <- sapply(1:num.markers, function(m){
    from <- cum.allele.counts[m]+1
    to <- cum.allele.counts[m+1]
    f.marker <- freqs[from:to]
    # handle missing data
    a <- which.max(f.marker)
    f.marker[a] <- 1.0 - sum(f.marker[-a])
    sum(f.marker^2)
  })
  he <- mean(1 - f.square.sums)
  return(he)
}

# E-E, E-NE and A-NE evaluation of core (genotypes, phenotypes or precomputed distances)
# Uses Modified Rogers for genotypes and Gower for phenotypes
ee <- function(core, data){
  meas <- get.measure(data)
  evaluateCore(core, data, objective("EE", meas))
}
en <- function(core, data){
  meas <- get.measure(data)
  evaluateCore(core, data, objective("EN", meas))
}
an <- function(core, data){
  meas <- get.measure(data)
  evaluateCore(core, data, objective("AN", meas))
}
get.measure <- function(data){
  if(is(data, "chgeno")){
    meas <- "MR"
  } else if(is(data, "chpheno")){
    meas <- "GD"
  } else if(is(data, "chdist")){
    meas <- "PD"
  } else {
    stop("Data should be genotypes, phenotypes or a precomputed distance matrix.")
  }
  return(meas)
}

# minimum distance evaluation (precomputed distance matrix)
dmin <- function(core, dist){
  if(is(core, "chcore")){
    core <- core$sel
  }
  if(is(dist, "chdist")){
    dist <- dist$data
  }
  dist.core <- dist[core, core]
  min.dist <- min(dist.core[upper.tri(dist.core)])
  return(min.dist)
}

# allelic diversity evaluations (genotypes only)
he <- function(core, geno){
  evaluateCore(core, geno, objective("HE"))
}
sh <- function(core, geno){
  evaluateCore(core, geno, objective("SH"))
}