#' Turn wide/matrix RFI data into tidy data
#' 
#' Takes wide/matrix RFI data exported from Zeptosens where each row is an
#' sample and each column is an antibody and turns it into \href{https://cran.r-project.org/web/packages/tidyr/vignettes/tidy-data.html}{tidy data}.
#' 
#' A "Batch" column is added by default so that batches can be identified if
#' tidy dataframes from several runs are joined together. 
#' 
#' If the exact same sample name appears multiple times, they represent 
#' technical replicates and can be averaged using 'ave_reps'.
#' 
#' 
#' @return Tidy dataframe with columns containing full antibody names and 
#' phenotype information can also be merged into the base RFI dataframe. 
#' 
#' 
#' @param df Dataframe in the (matrix) format of the raw data.
#' @param ABnames Dataframe containing the antibody number - full
#'     antibody name key-pairs.
#' @param Batch Batch code of the run. Can be single string or number. A column
#'     of this code will be added to the output dataframe.
#' @param ave_reps Single logical indicating whether technical replicates 
#'     (samples with the same names) should be averaged.
#' @param pheno Dataframe containing information on the sample phenotypes.
#' 
#' @importFrom assertthat assert_that
#' @importFrom dplyr mutate
#' @importFrom dplyr summarise
#' @importFrom dplyr group_by
#' 
#' 
#' @export
tidyData <- function(df, ABnames, Batch = "A", ave_reps = FALSE, pheno){
  
  # check argument inputs
  assert_that(colnames(df)[1] == "X1", dim(df)[2] > 1,
                          msg = "Check 'df' dataframe")
  
  assert_that(
    sum(c("Antibody.Name","Ab.No.") %in% colnames(ABnames)) == 2,
    dim(ABnames)[2] == 2,
    msg = "Check column names in 'ABnames' dataframe")
  
  assert_that(length(Batch) == 1,
                          msg = "Batch should have length of 1")
  
  assert_that(is.logical(ave_reps),
                          length(ave_reps) == 1,
                          msg = "Check 'ave_reps' is a single logical")
  
  if (! missing(pheno)){
    assertthat::assert_that(sum(colnames(pheno) %in% c("Lysate.ID")) == 1,
                            msg = "Check pheno dataframe has one 'Lysate.ID' column")
  }
  
  
  # gather data
  numcols <- ncol(df)
  
  gather_df <- df %>%
    tidyr::gather(2:numcols, key = "AB", value = "RFI") %>%
    mutate(Batch = Batch)
  
  # average technical replicates
  if (ave_reps){
    
    gather_df <- gather_df %>%
      group_by(AB,X1) %>%
      summarise(RFI = mean(RFI, na.rm = TRUE)) %>%
      mutate(Batch = Batch)
    
  } 
  
  
  # merge ABnames data
  tidydf <- merge(gather_df, ABnames, by.x = "AB", by.y = "Ab.No.",
                  all.x = TRUE, all.y = FALSE)
  
  # merge pheno data, if input given
  if (! missing(pheno)){
    
    tidydf <- merge(tidydf, pheno, by.x = "X1", by.y = "Lysate.ID",
                    all.x = TRUE, all.y = FALSE)
    
  }
  
  return(tidydf)
  
}