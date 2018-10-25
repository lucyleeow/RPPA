#' Convert tidy dataframe into matrix
#' 
#' Converts a tidy dataframe into a numeric matrix.
#' 
#' If the tidy dataframe contains technical replicates with the same sample
#' names, unique sample names will be created by add a group row number to the
#' the end of the sample name.
#' 
#' Useful for functions that require data to be in matrix format like 
#' \code{pheatmap} and \code{plotMA}.
#' 
#' @return A numeric matrix where each row is an antibody and each column is
#'     a sample.
#' 
#' @param tidydf Tidy dataframe to convert to matrix 
#' @param logdata Single logical indicating whether the raw RFI values should
#'     be logged.
#' @param tech_reps Single logical indicating whether ther are technical
#'     replicates.
#'     
#' @importFrom assertthat assert_that
#' 
#' 
#' @export
df_to_mat <- function(tidydf, logdata, tech_reps){
  
  # check inputs
  assert_that(is.logical(logdata),
                          length(logdata) == 1,
                          msg = "Check 'logdata' is a single logical")
  
  assert_that(is.logical(tech_reps),
                          length(tech_reps) == 1,
                          msg = "Check 'tech_reps' is a single logical")
  
  
  # create unqiue sample names if there are technical replicates
  if (tech_reps){
    
    tidydf <- tidydf %>%
      group_by(X1, AB) %>%
      mutate(X2 = paste(X1, row_number(), sep = "_")) %>%
      ungroup() %>%
      select(-X1) %>%
      rename(X1 = X2)
    
  }
  

  # make matrix
  mat <- as.matrix(
    tidydf %>%
      select(X1, AB, RFI) %>%
      spread(key = AB, value = RFI)
  )
  
  rownames(mat) <- mat[,1]
  mat <- mat[,-1]
  mode(mat) <- "numeric"
  mat <- t(mat)
  
  # log result adding prior of 0.00001 to avoid taking log of 0
  
  if (logdata){
    mat <- log2(mat + 0.00001)
  }
  
  return(mat)
  
}