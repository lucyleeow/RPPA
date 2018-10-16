###################################################
# Functions to clean/format RPPA data from excel
# Author: Lucy Liu
###################################################

############################
# Comparison 
############################


############################
# tidyIinput 
############################ 
#
# This function takes input RPPA data of format:
# Col1          Col2..... Colx
# Sample name   AB1       ABn
# 
# and:
# 1. adds antibody name
# 2. converts format into tidy data.
#

tidyInput <- function(
  df,           # name of the df to turn into tidy format
  ABnames,      # df containing the AB code-name conversions
  Batch = "A",  # Batch code of the run
  reps = FALSE, # logical indicating whether there were technical. Replicates which will be averaged
  pheno         # optional argument. Name of the phenotype df to be merged      
){
  
  # check argument inputs
  
  ## check df  
  ## check that the first column is called 'X1' and that there are >1 cols
  assertthat::assert_that(colnames(df)[1] == "X1", dim(df)[2] > 1,
                          msg = "Check 'df' dataframe")
  
  ## check ABnames
  assertthat::assert_that(
    sum(colnames(ABnames) %in% c("Antibody.Name","Ab.No.")) == 2,
    dim(ABnames)[2] == 2,
    msg = "Check column names in 'ABnames' dataframe")
  
  ## check ABnames
  assertthat::assert_that(length(Batch) == 1,
                          msg = "Batch should be a vector of 1")
  
  ## check reps input
  assertthat::assert_that(is.logical(reps), length(reps) == 1,
              msg = "Check 'reps' is a single logical.")

  ## check pheno
  assertthat::assert_that(sum(colnames(pheno) %in% c("Lysate.ID")) == 1,
                              msg = "Check pheno dataframe has one 'Lysate.ID' column")

  # gather data
  
  numcols <- ncol(df)
  
  gather_df <- df %>%
    gather(2:numcols, key = "AB", value = "RFI") %>%
    mutate(Batch = Batch)
  
  if (reps){
    
    gather_df <- gather_df %>%
      group_by(AB,X1) %>%
      summarise(RFI = mean(RFI, na.rm = TRUE)) %>%
      mutate(Batch = Batch)
    
  } 
  
  # merge ABnames data
  tidydf <- merge(gather_df, ABnames, by.x = "AB", by.y = "Ab.No.")
  
  # merge pheno data, if input given
  if (! missing(pheno)){
    
    tidydf <- merge(tidydf, pheno, by.x = "X1", by.y = "Lysate.ID")
    
  }
  
  return(tidydf)
  
}

############################
# matIinput 
############################# 
#
# This function turns tidy dataframe of RPPA data into a matrix, 
# with row = AB's and columns = Samples. Rownames are AB codes
# 

matInput <- function(
  tidydf,       # name of the tidydf to turn into matrix format
  logdata      # logical indicating whether the RFI values should be log2'ed
){
  
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

