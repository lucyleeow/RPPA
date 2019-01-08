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
#' @param tidydf Tidy dataframe of RFI values to convert to matrix.
#' @param logdata Single logical indicating whether the raw RFI values should
#'     be logged.
#' @param tech_reps Single logical indicating whether there are technical
#'     replicates.
#' @param convert_neg Single logical indicating whether negative values should
#'     be converted to zero (0).
#' @param AB_col Name of the column containing antibody name or number,
#'     as string. These names or numbers will become rownames of the 
#'     output matrix.
#'    
#'     
#' @importFrom assertthat assert_that
#' @importFrom magrittr %>%
#' 
#' 
#' @export
df_to_mat <- function(tidydf, logdata, tech_reps, convert_neg = TRUE,
                      AB_col = "Antibody.Name") {
  
  # check inputs
  assert_that(is.logical(logdata),
                          length(logdata) == 1,
                          msg = "Check 'logdata' is a single logical")
  
  assert_that(is.logical(tech_reps),
                          length(tech_reps) == 1,
                          msg = "Check 'tech_reps' is a single logical")
  
  assert_that(is.logical(convert_neg),
              length(convert_neg) == 1,
              msg = "Check 'convert_neg' is a single logical")
  
  
  # create unique sample names if there are technical replicates
  if (tech_reps){
    
    tidydf <- tidydf %>%
      dplyr::group_by(X1, !!as.name(AB_col)) %>%
      dplyr::mutate(X2 = paste(X1, row_number(), sep = "_")) %>%
      dplyr::ungroup() %>%
      dplyr::select(-X1) %>%
      dplyr::rename(X1 = X2)
    
  }
  

  # make matrix
  mat <- as.matrix(
    tidydf %>%
      dplyr::select(!!as.name(AB_col), AB, RFI) %>%
      tidyr::spread(key = !!as.name(AB_col), value = RFI)
  )
  
  rownames(mat) <- mat[,1]
  mat <- mat[,-1]
  mode(mat) <- "numeric"
  mat <- t(mat)
  
  
  # convert negative values to 0
  if (convert_neg) {
    
    mat[mat<0] <- 0
    
  }
  
  # log result adding prior of 0.00001 to avoid taking log of 0
  
  if (logdata){
    mat <- log2(mat + 0.00001)
  }
  
  return(mat)
  
}