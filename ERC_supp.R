##supplementary ERC functions
##updated GH 10/28/24


##make sure the gene list is contained in the erc matrix
clean_list = function(list, names){
  list <- unique(list)
  list <- list[which(match(list,names) != "NA")]
  list
}

pair_list <- function(list, erc_matrix , na.val = -2) {
  #provide vector list of genes (row and col names) and erc_matrix
  #returns vector list of gene x gene ERC values with no repeats or 'self x self' entries
  # good for computing averages and making histograms, etc...
  list <- clean_list(list , colnames(erc_matrix))
  mat <- erc_matrix[list,list]
  mat[which(mat == na.val)] <- NA
  mat[lower.tri(mat)]
}
##Make the ERC matrix symmetrical
make_symmetric = function(erc_matrix){
  size = dim(erc_matrix)[1]
  hold = matrix(NA,size,size)
  nindex = is.na(erc_matrix)   # find NA indices, make them zero
  hold[nindex] <- 0
  hold <- hold + t(erc_matrix) # add transpose to zeroes
  erc_matrix[nindex] <- hold[nindex] # combine matrices
  erc_matrix
}

##get permutations for a list of genes against the genome
permTestMat <- function(list, erc_matrix, perms=10000) {
  cleannames <- clean_list(colnames(erc_matrix),colnames(erc_matrix))
  list <- clean_list(list, cleannames)
  if(length(list) == 1){return(NaN)}
  obs <- mean(erc_matrix[list,list] , na.rm=TRUE)
  print(obs)
  if (is.nan(obs)) { return(NaN) }
  
  null <- c()
  while (length(null) < perms) {
    s <- sample(cleannames , length(list))
    m <- mean(pair_list(s , erc_matrix) , na.rm=TRUE)
    if (is.nan(m)) next
    null <- c(null , m)
  }
  pval <- sum(as.numeric(obs <= null)) / perms
  list( obs=obs , p=pval , null=null)
}
##Supply the output from computeERC to get the Fisher transformed correlation values
fishertransformed_updated <- function(erc_matrix){
  out <- erc_matrix$cor
  out <- atanh(erc_matrix$cor)*sqrt((erc_matrix$count)-3)
  ##out[upper.tri(out)] <- erc_matrix[upper.tri(erc_matrix)]
  out[is.nan(out)] <- NA
  out
}

##Get a subset matrix based on two gene lists
betweencomplex <- function(list1 , list2 , erc_matrix) {
  list1 <- clean_list(list1 , colnames(erc_matrix))
  list2 <- clean_list(list2 , colnames(erc_matrix))
  mat <- erc_matrix[list1,list2]
  mat
}
