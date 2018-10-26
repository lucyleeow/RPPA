#' Pairwise correlation coefficient plot
#' 
#' Calculate pairwise Spearman correlation between all possible combinations
#' of proteins and plot these (ordered) as a scatter plot.
#' 
#' @param 
#' 




################################################
# Helper function - gives vector of correlations
################################################
rankGR <- tidyGR %>%
  group_by(AB) %>%
  mutate(Rank = rank(RFI)) %>%
  select(-RFI)


spearmanRanks <- function(
  rankdf,    # Dataframe containing column of AB ranks and names of ABs
  c1         # Name of the column containg the ranks
  )
  {
  
  # Takes in a df where there is a column called "AB" 
  # (which contains the names of ABs) and another column of ranks of 
  # each AB (across all samples) and returns a vector of spearman correlations 
  # between the ranks of all possible combinations of pairs of ABs.
  rankdf <- tidydf %>%
    group_by(AB) %>%
    mutate(Rank = rank(RFI)) %>%
    select(-RFI)
  

    ABcomb <- combn(unique(rankdf$AB),2)
  # Performs "(Number of unique AB) Choose 2" and gives all the possible
  # combinations in matrix format. Matrix will have 2 rows as we are 
  # and "(Number of unique AB) Choose 2" columns.
  
  c <- vector(mode = "double", length = dim(ABcomb)[2])
  # Create empty vector length of the number of combinations
  
  for (i in 1:dim(ABcomb)[2]){
    # Cycle i through the number of combinations
    
    a <- rankdf[[c1]][rankdf$AB == ABcomb[1,i]]
    # Note cannot use rankdf$c1 as it will look for the column named "c1". Used instead rankdf[[c1]]. Double brackets returns vector instead of df/tibble
    b <- rankdf[[c1]][rankdf$AB == ABcomb[2,i]]
    
    c[i] <- cor(a,b, method = 'spearman')
  }
  
  return(c)
  
}

################################
# Main function - produces plot
################################

spearmanPlot <- function(rankdf, c1, legend){
  
  a <- spearmanRanks(rankdf = rankdf, c1 = as.character(c1))
  # gives you vector of correlations for the normalisation method of interest
  b <- spearmanRanks(rankdf = rankdf, c1 = "RFI")
  # gives you vector of correlations for the raw RFI values
  plot(a[order(a)], ylim = c(-1,1))
  points(b[order(b)], col = "Blue")
  legend("bottomright", legend = c(legend, "Raw"), text.col = c("Black", "Blue"))
  abline(0.5,0)
  abline(-0.5,0)
  # plot the ordered values.
  
}