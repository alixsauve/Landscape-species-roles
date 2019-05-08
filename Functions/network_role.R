# membership() returns a matrix b of belonging coefficients, dim(b) = S*H
# bij is the belonging coefficient of species i to community j.
membership <- function(t){
  # t is a 3d-array describing adjacencies in the landscape
  # dim(t) = S*S*H, S=number of species, H=number of habitats
  if (length(dim(t)) != 3){ # if t is not a tensor.
    stop("t must be a tensor.")
  }
  if (dim(t)[1] != dim(t)[2]){ # if each layer of t is not a tensor.
    stop("t must be an adjacency tensor.")
  }
  S <- dim(t)[1] # species diversity
  H <- dim(t)[3] # number of habitats
  b <- matrix(0, S, H) # matrix of belonging coefficients
  wh <- apply(t, c(1, 3), sum) # species weight in each habitat, dim(w) = S*H
  w <- apply(t, 1, sum) # species total weight
  b <- wh # initialising b
  for (i in 1:H){
    b[, i] <- b[, i]/w
  }
      
  return(b)
}

# c_value() returns the participation coefficient of a species to modules = among-module connectivity
# This function calls membership().
c_value <- function(t){
  # t is a 3d-array describing adjacencies in the landscape
  # dim(t) = S*S*H, S=number of species, H=number of habitats
  if (length(dim(t)) != 3){ # if t is not a tensor.
    stop("t must be a tensor.")
  }
  if (dim(t)[1] != dim(t)[2]){ # if each layer of t is not a tensor.
    stop("t must be an adjacency tensor.")
  }
  S <- dim(t)[1] # number of species
  H <- dim(t)[3] # number of habitats
  C <- rep(1, S) # partition coefficient of all species // initialisation
  beta <- membership(t) # belonging coefficient of species to each habitat
  for (i in 1:H){
    C <- C-beta[, i]^2
  }
  return(C)
}


# z_value() returns the within-module degree/weight of a species to a given module
z_value <- function(t, n){
  # t is a 3d-array describing adjacencies in the landscape
  # dim(t) = S*S*H, S=number of species, H=number of habitats
  # n is the module considered
  if (length(dim(t)) != 3){ # if t is not a tensor.
    stop("t must be a tensor.")
  }
  if (dim(t)[1] != dim(t)[2]){ # if each layer of t is not a tensor.
    stop("t must be an adjacency tensor.")
  }
  S <- dim(t)[1]
  H <- dim(t)[3]
  if ((n < 1) | (n > H)){
    stop("n is not a module/habitat number.")
  }
  Z <- vector("numeric", S) # within-module degree/weight of species
  W <- rowSums(t[,, n]) # species degree/weight in habitat n
  if (sd(W[W != 0]) != 0){
    Z[W != 0] <- scale(W[W != 0], center = TRUE, scale = TRUE)
  }
  Z[W == 0] <- NA
  return(Z)
}

# z_value_landscape() returns returns the mean within-habitat degree/weight of all species at the landscape scale
# this functions calls the following function: z_value()
z_value_landscape <- function(t){
  # t is a 3d-array describing adjacencies (binary or quantitative) in the landscape
  # dim(t) = S*S*H, S=number of species, H=number of habitats
  if (length(dim(t)) != 3){ # if t is not a tensor.
    stop("t must be a tensor.")
  }
  if (dim(t)[1] != dim(t)[2]){ # if each layer of t is not a tensor.
    stop("t must be a tensor describing adjacencies.")
  }
  S <- dim(t)[1] # number of species
  H <- dim(t)[3] # number of habitats
  Z <- vector("numeric", S) # z values vector // initialisation
  beta <- membership(t) # belonging coefficients
  for (i in 1:H){ # for each habitat
    sp_i <- which(rowSums(t[,, i]) != 0) # list (numbers) of species present in habitat i
    z_i <- z_value(t, i) # within-habitat degree/weight
    Z[sp_i] <- Z[sp_i] + beta[sp_i, i]*z_i[sp_i]
  }
  return(Z)
}

# classification() returns the landscape role of a species given its z-c coordinates
classification <- function(z, c, z_lim, c_lim){
  # z is the within-habitat degree/weight
  # c is the among-habitat connectivity
  # z_lim is the threshold for z distinguishing peripheral/connectors from habitat/landscape hubs
  # c_lim is the threshold for c distinguishing peripheral/habitat hubs from connectors/landscape hubs
  role <- NA
  if ((z < z_lim) & (c <= c_lim) & (c >= 0)){
    role <- "peripheral"
  }
  if ((z < z_lim) & (c > c_lim) & (c <= 1)){
    role <- "connector"
  }
  if ((z >= z_lim) & (c <= c_lim) & (c >= 0)){
    role <- "habitat_hub"
  }
  if ((z >= z_lim) & (c > c_lim) & (c <= 1)){
    role <- "landscape_hub"
  }
  return(role)
}
