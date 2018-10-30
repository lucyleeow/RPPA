#' Pairwise correlation coefficient plot
#' 
#' Calculate the pairwise Spearman correlation between all possible combinations
#' of proteins and plot these (ordered) as a scatter plot.
#' 
#' This function performs the following steps:
#' 
#' \enumerate{
#'     \item For each protein, get the rank of it's expression level across
#'         all samples.
#'     \item Get matrix of all possible pairs of proteins.
#'     \item For each pair of proteins, calculate the Spearman correlation
#'         of their ranks.
#'     \item Plot all the (ordered) Spearman correlations for all the protein
#'         pairs as a scatter plot.
#' }
#' 
#' The distribution of the correlations can be indicative of the performance of
#' a normalisation technique. Differences in sample loading can result in 
#' artifically high correlations between protein pairs. You expect well 
#' normalised data to have similar positive and negative correlation 
#' coefficients. More information can be found in 
#' \href{http://journals.sagepub.com/doi/pdf/10.4137/CIN.S13329}{Liu et al.}.
#' 
#' 
#' @param tidydf Tidy dataframe with RFI values, antibody names, sample names
#'     and any other optional columns.
#' @param RFI1 Name of the column of the RFI values to rank, as string.
#' @param RFI2 Optional argument. Name of the second column of RFI values to
#'     rank, as string.
#' 
#' 
#' @importFrom assertthat assert_that
#' @importFrom magrittr %>%
#' @importFrom rlang !!
#' 
#' @export
spearmanPlot <- function(tidydf, RFI1, RFI2){
  
  # check inputs
  assert_that(is.character(RFI1), length(RFI1) == 1,
              msg = "Check 'RFI1' is a single string")
  
  assert_that(is.character(RFI2), length(RFI2) == 1,
              msg = "Check 'RFI2' is a single string")
  
  if (! missing(RFI1)){
    assert_that(sum(c("AB", RFI1, RFI2) %in% colnames(tidydf)) == 3,
                msg = "Check 'tidydf' contains the columns 'AB' and your input 'RFI1' & 'RFI2'")
  } else {
    assert_that(sum(c("AB", RFI1) %in% colnames(tidydf)) == 2,
                msg = "Check 'tidydf' contains the columns 'AB' and your input 'RFI1'")
  }
  
  # generate ranks
  rankdf <- tidydf %>%
    group_by(AB) %>%  # each group contains RFI values for 1 AB but across 
                      # ALL samples. Thus the rank is across samples.
    mutate(Rank1 = rank(!!(as.name(RFI1))))
  
  if (! missing(RFI2)){
    
    rankdf <- rankdf %>%
      group_by(AB) %>%
      mutate(Rank2 = rank(!!(as.name(RFI2))))
    
  }
  
  # generate vector of correlations
  corr1 <- .spearmanRanks(rankdf, "Rank1")
  
  if (! missing(RFI2)) {
    
    corr2 <- .spearmanRanks(rankdf, "Rank2")
    
  }
  
  plot(corr1[order(corr1)], ylim = c(-1,1))
  
  if (! missing(RFI2)) {
    
    points(corr2[order(corr2)], col = "Blue")
    
  }
  
  # vector of legend names
  legend <- RFI1
  
  if (! missing(RFI2)){
    
    legend <- c(legend, RFI2)
    
  }
  
  legend("bottomright", legend = legend, text.col = c("Black", "Blue"))
  
  abline(0.5,0)
  abline(-0.5,0)
  
  
}

# Helper function to get vector of correlations between each possible pair of
# antibodies.

#' @keywords internal
.spearmanRanks <- function(rankdf, rank_col) {
  
  # get all possible combinations of AB's
  ABcomb <- combn(unique(rankdf$AB),2)
  # note combn() performs "unique(rankdf$AB) Choose 2" and gives all the 
  # possible combinations in matrix format. Our matrix will have 2 rows as we 
  # and "(Number of unique AB) Choose 2" columns.
  
  # total number of possible combinations
  num_combinations <- ncol(ABcomb)
  
  # create empty vector, the length of the number of combinations
  c <- vector(mode = "double", length = num_combinations)
  
  ## fill vector with correlations
  for (i in 1:num_combinations){
    
    AB1_ranks <- rankdf[[rank_col]][rankdf$AB == ABcomb[1,i]]
    
    AB2_ranks <- rankdf[[rank_col]][rankdf$AB == ABcomb[2,i]]
    
    # note1: rankdf[[rank_col]] gives you a vector
    
    # note2: each vector of ranks are in the same order in terms of 
    # samples because the input data was tidy
    
    c[i] <- cor(AB1_ranks, AB2_ranks, method = 'spearman')
  }
  
  return(c)
  
}



  
  