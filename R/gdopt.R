
# Variant of GDOpt algorithm from "Statistical Techniques for
# Defining Reference Sets of Accesssions and Microsatellite Markers"
# (van Eeuwijk et al., 2011)
gdopt <- function(dist, size){
  
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
  
  # cluster with K-medoids algorithm
  clust <- pam(dist, diss = TRUE, k = size)
  
  # include medoids in core
  core <- clust$medoids
  return(core)
  
}
