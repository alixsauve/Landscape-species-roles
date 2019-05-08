# Shannon diversity of resources

shannon_div_diet <- function(m){
  # m is the adjacency matrix of the directed network
  
  if (ncol(m) != nrow(m)){
    stop("m must be a square matrix.")
  }
  if (isSymmetric(m)){
    stop("m must be the adjacency matrix of the directed network.")
  }
  
  S <- ncol(m) # number of species
  div_diet <- vector("numeric", S)
  basal_sp <- which(colSums(m) == 0) # basal species
  div_diet[basal_sp] <- NA
  
  for (sp in 1:S){
    if (!is.na(div_diet[sp])){
      diet <- which(m[, sp] != 0)
      p_sp <- m[diet, sp]/sum(m[diet, sp])
      for (k in 1:length(diet)){
        div_diet[sp] <- div_diet[sp]-p_sp[k]*log(p_sp[k])
        if (is.na(div_diet[sp])){
          stop()
        }
      }
    }
  }
  return(div_diet)
}