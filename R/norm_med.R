#' Median normalisation
#' 
#' Perform median normalisation on a tidy dataframe of RFI data.
#' 
#' This method assumes that the amount of all proteins measured in the 
#' experiment reflect the total protein amount of one sample. Thus, the 
#' median AB of a sample estimates total protein amount or 'sample loading'. 
#'
#' The normalisation calculation is: divide all raw values by the median value
#' of all AB signals of that sample is used to normalise the raw intensity
#' values for each AB. This may be biased when the number of ABs is <100.
#' 
#' @return The same input dataframe but with a 
#' 
#' @param tidydf Tidy dataframe of RFI values. There should be NO technical
#'     replicates.
#' @param col_name Column name of the median normalised values as a string.
#'     
#' @importFrom assertthat assert_that
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @importFrom dplyr select
#'
#' @export
norm_med <- function(tidydf, col_name = "RFI") {
  
  # check inputs
  assert_that(sum(c("X1", "AB", "RFI") %in% colnames(tidydf)) == 3,
                          msg = "Check 'tidydf' has 'X1', 'AB', 'RFI' columns")
  
  assert_that(nrow(unique(tidydf[,c("X1","AB")])) == nrow(tidydf),
                          msg = "Check that there are no technical replicates")
  
  assert_that(is.character(col_name), length(col_name) == 1,
              msg = "Check 'col_name' is a single string")
  
  
  # perform normalisation
  norm_df <- tidydf %>%
    group_by(X1) %>%
    mutate(med = median(RFI)) %>%
    dplyr::ungroup() %>%
    mutate(col_name = RFI/med) %>%
    select(-med)
  
  
  return(norm_df)
  
}