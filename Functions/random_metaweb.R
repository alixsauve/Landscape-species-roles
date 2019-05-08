# null models for random metawebs
my_sample <- function(set, n, repl, pi=NULL){
  # set is a vector of elements to sample
  # n is the number of elements to sample in set
  # repl is a logical, =TRUE if same element of set can be sampled several times
  # pi is a vector of probability to pick a given element from set
  if ((n > length(set)) & (repl == FALSE)){
    stop("Cannot sample more elements than present in set if repl = FALSE")
  }
  if (n<0){
    stop("n should be a postive integer.")
  }
  if ((length(set) != length(pi)) & (is.null(pi) == FALSE)){
    stop("pi and set should have the same length.")
  }
  res_sample <- numeric(0)
  
  if (length(set) != 1){
    res_sample <- sample(set, n, replace = repl, prob = pi)
  }
  else {
    res_sample <- rep(set, n)
  }
  return(res_sample)
}


random_mw <- function(mw, list_sp, lg_occ, X, lg_A = NULL){
  # mw is the 3d-array of adjacencies in the landscape, dim = S*S*H
  # list_sp is the list of species: guild, taxon, interaction mode
  # lg_occ is a binary matrix of plant species occurrences
  # X corresponds to interaction events per species
  # lg_A corresponds to lower guild abundance
  if (length(dim(mw)) != 3){
    stop("mw should be a 3d-array.")
  }
  if (dim(mw)[1] != dim(mw)[2]){
    stop("Each slice of mw should be a square matrix.")
  }
  if (is.null(lg_A) == TRUE){
    lg_A <- lg_occ
  }
  # make sure there is a single type of interaction
  S <- dim(mw)[1] # number of species
  H <- dim(mw)[3] # number of habitats
  
  m <- apply(mw, c(1, 2), sum) # aggregated network

  S_h <- vector("numeric", H) # diversity in each habitat
  S_h_null <- vector("numeric", H) # diversity in each habitat of the null metaweb
  
  L_h <- apply(mw, 3, sum)/2 # number of links in each habitat
  L_h_null <- vector("numeric", H) # number of links in each habitat of the null metaweb
  
  X_null <- vector("numeric", S) # interaction events per species in the random network
  X_null[is.na(X) == TRUE] <- NA
  
  mw_null <- array(0, dim = dim(mw))
  mw_pot <- array(0, dim = dim(mw)) # potential interactions in each habitat given occurrence of lower guild species
  # mw_pot is not symmetric
  hab_list_sp <- vector("list", H) # species list in each habitat
  hab_lg_list <- vector("list", H) # lower guild list in each habitat
  # enter lower guild species numbers
  for (i in 1:H){
    lg_h <- which(lg_occ[, i] == 1)
    hab_lg_list[[i]] <- lg_h
    mw_pot[lg_h,, i] <- m[lg_h, ]
  }
  m_null <- apply(mw_null, c(1, 2), sum)
  
  # connect all species with a given number of interaction events at least ONCE (insects)
  species <- which((X-X_null > 0) & (is.na(X) == FALSE) & (colSums(m_null) == 0))
  while(length(species) != 0){
    sp <- my_sample(species, 1, FALSE)
    if (list_sp$guild[sp] == "lower"){ # if species belongs to the lower guild
      h_set <- which(lg_occ[sp, ] != 0)
      pred <- which((m[sp, ] != 0) & ((X-X_null) > 0)) # predators which can interact
      int_set <- expand.grid(pred, h_set)
      int <- sample(1:dim(int_set)[1], 1, FALSE) # pick one interaction
      mw_null[sp, int_set[int, 1], int_set[int, 2]] <- 1
      mw_null[int_set[int, 1], sp, int_set[int, 2]] <- 1
      X_null[sp] <- X_null[sp]+1
      X_null[int_set[int, 1]] <- X_null[int_set[int, 1]]+1
      
      m_null <- apply(mw_null, c(1, 2), sum)
      species <- which((X-X_null > 0) & (is.na(X) == FALSE) & (colSums(m_null) == 0))
    }
    else { # if species belongs to the upper guild
      if (unique(is.na(X[list_sp$guild == "lower"])) == FALSE){ # if lower taxa have a given number of interaction events
        diet <- which((m[sp, ] != 0) & ((X-X_null) > 0))
      }
      else {
        diet <- which(m[sp, ] != 0)
      }
      prey <- my_sample(diet, 1, FALSE)
      h_set <- which(lg_occ[prey, ] != 0)
      h <- my_sample(h_set, 1, FALSE, lg_A[prey, h_set])
      mw_null[sp, prey, h] <- 1
      mw_null[prey, sp, h] <- 1
      X_null[sp] <- X_null[sp]+1
      
      m_null <- apply(mw_null, c(1, 2), sum)
      species <- which((X-X_null > 0) & (is.na(X) == FALSE) & (colSums(m_null) == 0))
      if (is.na(X[prey]) == FALSE){
	      X_null[prey] <- X_null[prey]+1
      }
    }
  }

  m_null <- apply(mw_null, c(1, 2), sum)
  # some species may remain to be connected
  isolated <- which(colSums(m_null) == 0)
  while (length(isolated) > 0){
    sp <- my_sample(isolated, 1, FALSE) # pick an isolated species
    isolated <- isolated[-which(isolated == sp)]
    pred_set <- which((m[sp, ] != 0) & ((X-X_null) > 0)) # set of predators feeding on sp
    pred <- my_sample(pred_set, 1, FALSE) # pick a predator
    h_set <- which(lg_occ[sp, ] != 0)
    h <- my_sample(h_set, 1, FALSE)
    mw_null[sp, pred, h] <- 1
    mw_null[pred, sp, h] <- 1
    X_null[pred] <- X_null[pred]+1
    X_null[sp] <- X[sp]+1
  }
    
  m_null <- apply(mw_null, c(1, 2), sum)
  m_null[m_null != 0] <- 1 # qualitative random network
  m[m != 0] <- 1 # qualitative observed network
  # have all possible interactions been created once?
  d <- m-m_null
  int2create <- which((d != 0) & (upper.tri(d) == TRUE), arr.ind = TRUE) # non-zero elements of d correspond to interactions not created in mw_null
  L2create <- dim(int2create)[1] # number of interaction to create
  if (L2create != 0){
    for (i in 1:L2create){ # for each interaction to be created in mw_null
      lt_i <- int2create[i, 1] # lower taxon
      ut_i <- int2create[i, 2] # upper taxon
      lt_abund_ih <- lg_A[lt_i, ] # abundance of the lower guild in each habitat
      p_ih <- lt_abund_ih/sum(lt_abund_ih) # probability to create the interaction in each habitat (proportional to the lower taxon abundance)
      h_i <- my_sample(seq(1:H), 1, FALSE, p_ih)# select one habitat where the interaction is created
      mw_null[lt_i, ut_i, h_i] <- 1 # create interaction in h_i
      mw_null[ut_i, lt_i, h_i] <- 1
      X_null[c(lt_i, ut_i)] <- X_null[c(lt_i, ut_i)]+1 # add interaction event
    }
  }

  
  m_null <- apply(mw_null, c(1, 2), sum)

  # maybe some interaction events remain to be assigned
  # hereafter, we proceed by species from the upper guild to make sure we assign the right number of interactions per species
  species <- which(((X-X_null) > 0) & (list_sp$guild == "upper"))
  while (length(species) > 0){
    sp <- my_sample(species, 1, FALSE) # pick a species
    events <- (X-X_null)[sp] # number of interaction events to create
    mw_pot_sp <- array(0, dim = dim(mw))
    mw_pot_sp[sp,,] <- mw_pot[sp,,]
    mw_pot_sp[, sp,] <- mw_pot[, sp,]
    
    # weight potential interactions with lower guild abundance -> m_ijh = m_ijh*lg_A_ih
    lg_sp <- which(list_sp$guild == "lower")
    if (is.null(lg_A) == FALSE){
      for (i in 1:H){
        for (j in 1:S){
          mw_pot_sp[j, lg_sp, i] <- mw_pot_sp[j, lg_sp, i]*lg_A[, i]
          mw_pot_sp[lg_sp, j, i] <- mw_pot_sp[lg_sp, j, i]*lg_A[, i]
        }
      }
    }
    int_set <- which(mw_pot_sp != 0)
    int <- my_sample(int_set, events, TRUE, mw_pot_sp[int_set])
    freq <- summary(as.factor(int))
    int2 <- unique(int)
    mw_null[int2] <- mw_null[int2]+as.vector(freq)
    
    if (unique(is.na(X[list_sp$guild == "lower"])) == TRUE){
      X_null[list_sp$guild == "lower"] <- NA
    }
    X_null[sp] <- X_null[sp] + events
    species <- species[-which(species == sp)] # remove it from the set of species to connect more
  }

	# verify that all species in the upper guild have the right number of interaction events
  if (sum((X-X_null)[list_sp$guild == "upper"], na.rm=TRUE) != 0){
    print(X)
    print(X_null)
    stop()
  }
  
  for (i in 1:H){
    mw_i <- mw_null[,, i]
    if (isSymmetric(mw_i) == FALSE){
      mw_i[lower.tri(mw_i)] <- t(mw_i)[lower.tri(mw_i)]
      mw_null[,, i] <- mw_i
    }
  }
  m_null <- apply(mw_null, c(1, 2), sum)
  
  return(mw_null)
}

# full_random_mw() returns a random metaweb for all types of interactions using random_mw()
full_random_mw <- function(a, int, sp, occ, type_list, abund = NULL, f){
  # a is the observed quantitative metaweb with multiple types of interactions
  # int is the corresponding interaction list (mode, int_coord1, int_coord2)
  # sp is the species list
  # occ is the occurrence matrix, initially filled for plants
  # type_list is the list of interaction types present in the metaweb
  # abund is the vector of species abundances in each type of habitat
  # f is a character string for output format: = "int_list" or "3d-array"
  
  if (length(dim(a)) != 3){
    stop("a must be a 3d-array.")
  }
  if (length(which(a!=0))/2 != dim(int)[1]){
    stop("a and int must describe the same interactions.")
  }
  if (dim(occ)[1] != dim(sp)[1]){
    stop("occ must describe the occurrences of each species.")
  }
  if ((is.null(int$mode) == TRUE) | (is.null(int$int_coord1) == TRUE) | (is.null(int$int_coord2) == TRUE) | (is.null(int$lower_taxon) == TRUE)){
    stop("int should contain informations on each interaction such as: lower_taxon, mode, int_coord1, int_coord2.")
  }
  if ((is.null(sp$taxon) == TRUE) | (is.null(sp$guild) == TRUE)){
    stop("sp should contain informations on each species such as: guild (= {'plant', 'insect'}), taxon.")
  }
  if (length(type_list) != length(unique(int$mode))){
    stop("The metaweb should contain the same types of interactions than given in type_list.")
  }
  if ((f != "int_list") & (f != "3d-array")){
    stop("f is the format of the output: either 'int_list' or '3d-array'.")
  }
  
  int_null <- as.data.frame(matrix(0, 0, 7))
  names(int_null) <- c("habitat", "upper_guild", "lower_guild", "upper_taxon", "lower_taxon", "mode", "freq")
  habitats <- unique(int$habitat) # habitat list
  
  a_null <- array(0, dim = dim(a))
  for (i in 1:length(type_list)){
    a_i <- a
    int_i <- c(int$int_coord1[int$mode == type_list[i]], int$int_coord2[int$mode == type_list[i]]) # interactions of mode i
    a_i[-int_i] <- 0
    m_i <- apply(a_i, c(1, 2), sum)
    sp_i <- which(colSums(m_i) != 0) # species involved in interaction mode i
    a_i <- a_i[sp_i, sp_i, ] # corresponding sub-network
    m_i <- m_i[sp_i, sp_i]
    sp_list_i <- sp[sp_i, ]
    lg_names <- unique(int$lower_taxon[int$mode == type_list[i]])
    lg <- which(is.na(match(sp$taxon[sp_i], lg_names)) == FALSE)
    lg_occ_i <- occ[sp_i[lg], ]
    sp_list_i$guild[lg] <- "lower"
    sp_list_i$guild[which(sp_list_i$guild != "lower")] <- "upper"
    int_freq <- colSums(m_i)
    pl <- which(sp$guild[sp_i] == "plant")
    int_freq[pl] <- NA
    if (is.null(abund) == FALSE){
      lg_abund <- abund[sp_i[lg], , i]
    }
    else {
      lg_abund <- NULL
    }
    
    a_null_i <- random_mw(a_i, sp_list_i, lg_occ_i, int_freq, lg_abund) 
    occ_i <- apply(a_null_i, c(1, 3), sum)
    occ_i[occ_i != 0] <- 1
    occ[sp_i, ] <- occ_i
    a_null[sp_i, sp_i, ] <- a_null[sp_i, sp_i, ]+a_null_i
    
    if (sum(a_null_i) != sum(a_i)){
      stop("The number of interaction events is not right.")
    }

    if (f == "int_list"){
      ug_name <- unique(int$upper_guild[int$mode == type_list[i]])
      lg_name <- unique(int$lower_guild[int$mode == type_list[i]])
      int_null_i <- mw_list(a_null_i, sp_list_i, type_list[i], habitats, ug_name, lg_name)
      int_null <- rbind(int_null, int_null_i)
		if (sum(int_null_i$freq) != sum(int$freq[int$mode == type_list[i]])){
			stop()
		}
    }

    # calculate abundances of upper guild after creating random metaweb
    if (is.null(abund) == FALSE){
      ug_names <- sp_list_i$taxon[sp_list_i$guild == "upper"]
      ug <- which(is.na(match(sp$taxon, ug_names)) == FALSE)
      ug_i <- which(is.na(match(sp_list_i$taxon, ug_names)) == FALSE)
      for (j in 1:length(type_list)){
        abund[ug, , j] <- apply(a_null_i, c(1, 3), sum)[ug_i, ]
      }
    }    
  }
  
  if (f == "3d-array"){
    return(a_null)
  }
  else {
    return(int_null)
  }
}

