#' Create phenotype dataframe
#' 
#' Creates a dataframe containing phenotype information on the RPPA samples. 
#' 
#' @param path Single string containing the path to the phenotyle xlsx file.
#' @param sheet Single numeric indicating the sheet number or name of the 
#'     xlsx sheet to read in.
#' @param col_remove Optional argument. Vector of column number(s) to remove.
#' @param condition_cols Optional argument. Vector of column number(s) to
#'     paste together to create the 'Conditon' column.
#'     
#' @return Dataframe of phenotype information. A 'Condition' column of sample
#'     and treatment details concatenated together will be included if the 
#'     'condition_cols' argument is given. 
#'     

#' @importFrom assertthat assert_that


#' @export
inputPheno <- function(path, sheet, col_remove, condition_cols) {
  
  # check inputs
  assert_that(is.character(path), 
                          length(path) == 1,
                          msg = "Check that 'path' is single string")      
  
  if (! missing(col_remove)){
    assert_that(is.numeric(col_remove),
                            msg = "Check 'col_remove' is of integer type")
  }
  
  if (! missing(condition_cols)){
    assert_that(is.numeric(condition_cols),
                            msg = "Check 'condition_cols' is of integer type")
  }
  
  
  # read in sheet
  df <- openxlsx::read.xlsx(path, sheet = sheet)
  
  ## remove white spaces
  df <- data.frame(lapply(df, trimws, which = "both"))
  
  assert_that(sum(colnames(df) %in% c("Lysate.ID")) == 1,
                          msg = "Check that your sheet contains a 'Lysate.ID' column")
  
  
  # create condition column
  if (! missing(condition_cols)){
    
    df$Condition <- do.call(paste, c(df[,condition_cols],
                                     sep = "_"))
    
  }
  
  # remove unwanted columns
  if (! missing(col_remove)){
    
    df <- df[, -col_remove]
    
  }
  
  # remove duplicates
  df <- df[! duplicated(df$Lysate.ID),]
  
  
  return(df)
  
}


