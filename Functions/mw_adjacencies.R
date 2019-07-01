# mw_adjacencies() returns an array S*S*c (S = number of species, c = number of habitats)
mw_adjacencies <- function(intlist, sp, w, d, n){
  # intlist is the interaction list with the corresponding habitats
  # sp is the species list (taxon and guild)
  # w is logical variable: =TRUE if metaweb is weighted, =FALSE otherwise
  # d is logical variable: =TRUE if metaweb is oriented, =FALSE otherwise
  # n is logical variable: =TRUE if using normalised interaction frequencies, =FALSE otherwise
  if (!is.data.frame(intlist)){
    stop("intlist must be a data frame.")
  }
  if ((is.null(intlist$habitat)) | (is.null(intlist$lower_taxon)) | (is.null(intlist$upper_taxon)) |
      (is.null(intlist$freq))){
    stop("intlist must have at least 4 columns: 'habitat', 'lower_taxon', 'upper_taxon', and 'freq'.")
  }
  if ((!is.logical(w)) | (!is.logical(d)) | (!is.logical(n))){
    stop("w, d and n must be logical variables.")
  }
  if ((is.null(intlist$freq_norm)) & (n)){
    stop("If the array is to be filled with normalised interaction frequencies, intlist must have a column named 'freq_norm'.")
  }
  S <- dim(sp)[1]
  hab <- sort(unique(intlist$habitat))
  web <- array(0, c(S, S, length(hab)))
  
  for (i in 1:dim(intlist)[1]){
    habitat <- which(hab == intlist$habitat[i])
    sp1 <- which(sp$taxon == intlist$lower_taxon[i])
    sp2 <- which(sp$taxon == intlist$upper_taxon[i])
    if (w == FALSE){
      web[sp1, sp2, habitat] <- 1
      web[sp2, sp1, habitat] <- 1
    }
    else {
      if (n == TRUE){
        web[sp1, sp2, habitat] <- intlist$freq_norm[i]
        web[sp2, sp1, habitat] <- intlist$freq_norm[i]
      }
      else {
        web[sp1, sp2, habitat] <- intlist$freq[i]
        web[sp2, sp1, habitat] <- intlist$freq[i]
      }
    }
    if (d == TRUE){
      web[sp2, sp1, habitat] <- 0
    }
  }
  return(web)
}
