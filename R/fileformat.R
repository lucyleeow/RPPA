###################################################
# Functions to clean/format RPPA data from excel
# Author: Lucy Liu
###################################################

library(openxlsx)
library(tidyverse)


#####################################
# Create design matrix for comparisons
#####################################

inputPheno <- function(
  path,            # path of the phenotype xlsx file
  sheet,           # sheet of xlsx file to read in
  col_remove,      # optional argument. column number(s) to remove
  condition_cols   # optional argument. column number(s) to paste 
                   # together to create 'condition' column
){
  
  # check inputs
  assertthat::assert_that(is.character(path), 
                          length(path) == 1,
                          msg = "Check that 'path' is single string")      
  
  if (! missing(col_remove)){
    assertthat::assert_that(is.integer(col_remove),
                            msg = "Check 'col_remove' is of integer type")
  }
  
  if (! missing(condition_cols)){
    assertthat::assert_that(is.integer(condition_cols),
                            msg = "Check 'condition_cols' is of integer type")
  }
  
  
  # read in sheet
  df <- read.xlsx(path, sheet = sheet)
  
  assertthat::assert_that(sum(colnames(df) %in% c("Lysate.ID")) == 1,
                          msg = "Check that your sheet contains a 'Lysate.ID' column")
  
  
  # remove unwanted columns
  if (! missing(col_remove)){
    
    df <- df[, -col_remove]
    
  }
  
  
  # create condition column
  if (! missing(condition_cols)){
    
    df$condition <- do.call(paste, c(df[,condition_cols],
                                     sep = "_"))
    
  }

  # remove duplicates
  df <- df[! duplicated(df$Lysate.ID),]
  
  
  return(df)
  
}


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

tidyData <- function(
  df,               # name of the df to turn into tidy format
  ABnames,          # df containing the AB code-name conversions
  Batch = "A",      # Batch code of the run
  ave_reps = FALSE, # logical indicating whether technical 
                    # replicates should be averaged
  pheno             # optional argument. Name of the phenotype df to be merged      
){
  
  # check argument inputs
  
  assertthat::assert_that(colnames(df)[1] == "X1", dim(df)[2] > 1,
                          msg = "Check 'df' dataframe")
  
  assertthat::assert_that(
    sum(colnames(ABnames) %in% c("Antibody.Name","Ab.No.")) == 2,
    dim(ABnames)[2] == 2,
    msg = "Check column names in 'ABnames' dataframe")

  assertthat::assert_that(length(Batch) == 1,
                          msg = "Batch should have length of 1")

  assertthat::assert_that(is.logical(ave_reps),
                          length(ave_reps) == 1,
                          msg = "Check 'ave_reps' is a single logical")

  if (! missing(pheno)){
    assertthat::assert_that(sum(colnames(pheno) %in% c("Lysate.ID")) == 1,
                              msg = "Check pheno dataframe has one 'Lysate.ID' column")
  }
  

  # gather data
  numcols <- ncol(df)
  
  gather_df <- df %>%
    gather(2:numcols, key = "AB", value = "RFI") %>%
    mutate(Batch = Batch)
  
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

############################
# matInput 
############################# 
#
# This function turns tidy dataframe of RPPA data into a matrix, 
# with row = AB's and columns = Samples. Rownames are AB codes
# 

matInput <- function(
  tidydf,       # name of the tidydf to turn into matrix format
  logdata,      # logical indicating whether the RFI values should be log2'ed
  tech_reps     # logical indicating whether there are technical
                # replicates (which are named the same)
){
  
  # check inputs
  assertthat::assert_that(is.logical(logdata),
                          length(logdata) == 1,
                          msg = "Check 'logdata' is a single logical")
  
  assertthat::assert_that(is.logical(tech_reps),
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

