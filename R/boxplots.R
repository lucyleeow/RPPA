#' Generate boxplots of RFI values
#' 
#' Generate boxplots of raw RFI value for each antibody or each sample. 
#' 
#' @param tidydf Tidy dataframe of raw RFI values.
#' @param RFIcol Name of the column, as string, containing RFI values.
#' @param fill_by Optional argument. Name of the column, a string, to add to
#'     the \code{fill} argument in ggplot's \code{aes}.
#'     
#'     
#' @importFrom assertthat assert_that
#' 
#' @describeIn RFIperAB Generates a boxplot for each antibody.
#' @export
RFIperAB <- function(tidydf, RFIcol = "RFI", fill_by){
  
  assert_that(RFIcol %in% colnames(tidydf) == 1,
              msg = "Check 'RFIcol' is correct")
  
  if (! missing(fill_by)){
    assert_that(fill_by %in% colnames(tidydf) == 1,
                msg = "Check 'fill_by' is correct")
  }
  
  assert_that("Antibody.Name" %in% colnames(tidydf) == 1,
              msg = "Check that an 'Antibody.Name' column exists in tidydf")
  
  
  # make plot
  if (missing(fill_by)){
    
    gg <- ggplot(tidydf, aes_string(y=RFIcol, x="Antibody.Name"))
    
  } else{
    
    gg <- ggplot(tidydf, aes_string(y=RFIcol, x="Antibody.Name", fill=fill_by))
  }
  
  # make rest of plot
  
  gg +
    geom_boxplot() + 
    labs(title = "RFI per AB", y = "RFI", x = "Antibody") +
    coord_flip() +
    theme(plot.title = element_text(hjust = 0.5), 
          title = element_text(size=16),
          axis.text.y = element_text(size = 11))
  
}

#' @describeIn RFIperAB Generates a boxplot for each sample.
#' @export
RFIperSample <- function(tidydf, RFIcol = "RFI", fill_by){
  
  # check inputs
  assert_that(RFIcol %in% colnames(tidydf) == 1,
              msg = "Check 'RFIcol' is correct")
  
  if (! missing(fill_by)){
    assert_that(fill_by %in% colnames(tidydf) == 1,
                msg = "Check 'fill_by' is correct")
  }
  
  assert_that("Antibody.Name" %in% colnames(tidydf) == 1,
              msg = "Check that an 'Antibody.Name' column exists in tidydf")
  
  
  if (missing(fill_by)){
    
    gg <- ggplot(tidydf, aes_string(y=RFIcol, x="X1"))
    
  } else{
    
    gg <- ggplot(tidydf, aes_string(y=RFIcol, x="X1", fill=fill_by))
  }
  
  # make plot
  
  gg + 
    geom_boxplot() + 
    labs(title = "RFI per Sample", y = "RFI", x = "Sample") +
    coord_flip() +
    theme(plot.title = element_text(hjust = 0.5), 
          title = element_text(size=16),
          axis.text.y = element_text(size = 11))
}