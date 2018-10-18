#################################################################
# Functions to normalise RPPA data using various normalisation
# methods
#
# Author: Lucy Liu
#################################################################


##########################################################
# Perform median normalisation
#
# This method assumes that the amount of all proteins measured in the 
# experiment reflect the total protein amount of one sample. Thus, the 
# median AB of a sample estimates sample loading. 
#
# Divide all raw values by the median value of all AB signals of that
# sample is used to normalise the raw intensity values for each AB. 
# May be biased when the number of ABs is <100
##########################################################


norm_med <- function(
  df    # df in tidy format with NO technical replicates
){
  
  # check inputs
  assertthat::assert_that(sum(c("X1", "AB", "RFI") %in% colnames(df)) == 3,
                              msg = "Check 'df' has 'X1', 'AB', 'RFI' columns")
  
  assertthat::assert_that(nrow(unique(df[,c("X1","AB")])) == nrow(df),
                          msg = "Check that there are no technical replicates")
  
  
  norm_df <- df %>%
    group_by(X1) %>%
    mutate(med = median(RFI)) %>%
    ungroup() %>%
    mutate(RFI = RFI/med) %>%
    select(-med)
  
  
  return(norm_df)
  
}
