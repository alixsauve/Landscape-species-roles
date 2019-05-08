# mw_adjacencies() returns an array S*S*c (S = number of species, c = number of habitats)
mw_adjacencies <- function(list, sp, w, d, n){
  # list is the interaction list with the corresponding habitats
  # sp is the species list (taxon and guild)
  # w is logical: =TRUE if metaweb is weighted, =FALSE otherwise
  # d is logical: =TRUE if metaweb is oriented, =FALSE otherwise
  # n is logica: =TRUE if using normalised interaction frequencies, =FALSE otherwise
  
  if (is.logical(w) == FALSE){
    stop("w must be a logical.")
  }
  S <- dim(sp)[1]
  hab <- sort(unique(list$habitat))
  web <- array(0, c(S, S, length(hab)))
  
  for (i in 1:dim(list)[1]){
    habitat <- which(hab == list$habitat[i])
    sp1 <- which(sp$taxon == list$lower_taxon[i])
    sp2 <- which(sp$taxon == list$upper_taxon[i])
    if (w == FALSE){
      web[sp1, sp2, habitat] <- 1
      web[sp2, sp1, habitat] <- 1
    }
    else {
      if (n == TRUE){
        web[sp1, sp2, habitat] <- list$freq_norm[i]
        web[sp2, sp1, habitat] <- list$freq_norm[i]
      }
      else {
        web[sp1, sp2, habitat] <- list$freq[i]
        web[sp2, sp1, habitat] <- list$freq[i]
      }
    }
    if (d == TRUE){
      web[sp2, sp1, habitat] <- 0
    }
  }
  return(web)
}