#' Barplots per sample/antibody
#' 
#' Creates separate barplots of the RFI values for EACH sample or antibody and 
#' writes to either pdf or the current graphics device.
#' 
#' @param tidydf Tidy dataframe of RFI values to plot.
#' @param RFIcol Column name of the RFI column as string.
#' @param pdfoutput Single logical indicating whether the plots should be
#'     output to pdf file or the current graphics device.
#' 
#' 
#' @describeIn plotperSample Create a barplot of the RFI values of each 
#'     antibody. One barplot is created for each sample. Write to a file called 
#'     'plotPerSample.pdf' or the current graphics device.
#'     
#' @importFrom assertthat assert_that
#' @importFrom magrittr %>%
#' @import ggplot2
#' @importFrom dplyr filter 
#'     
#'  
#' @export
plotperSample <- function(tidydf, RFIcol = "RFI", pdfoutput = TRUE) {
  
  # check inputs
  assert_that(RFIcol %in% colnames(tidydf) == 1,
              msg = "Check 'RFIcol' exists in 'tidydf'")
  
  assert_that(is.logical(pdfoutput), length(pdfoutput) == 1,
              msg = "Check 'pdfoutput' is a single logical")
  
  
  # number of unique samples
  numSample <- length(unique(tidydf$X1))
  
  
  # make plot
  if (pdfoutput){
    pdf("RFI_barplots/plotPerSample.pdf")
  }
  
  for(i in 1:numSample){
    
    print(
      tidydf %>%
        filter(X1 == unique(tidydf$X1)[i]) %>%    # filter for only sample i
        ggplot(aes_string(y=RFIcol,x="Antibody.Name")) + 
        geom_bar(stat = "identity") +
        labs(title = unique(tidydf$X1)[i], y = "Average RFI", x = "Antibody") +
        theme(plot.title = element_text(hjust = 0.5),
              title = element_text(size=11), 
              axis.text.y = element_text(size = 6)) +
        coord_flip()
    )
  }
  
  if (pdfoutput){
    dev.off()
  }
  
}

#' @describeIn plotperSample Creates a barplot of the RFI values of each 
#'     sample. One barplot is created for each antibody. Writes to either a pdf called 
#'     'plotPerAB.pdf' or the current graphics device.
#' @export
plotperAB <- function(tidydf, RFIcol = "RFI", pdfoutput = TRUE) {
  
  # check inputs
  assert_that(RFIcol %in% colnames(tidydf) == 1,
              msg = "Check 'RFIcol' exists in 'tidydf'")
  
  assert_that(is.logical(pdfoutput), length(pdfoutput) == 1,
              msg = "Check 'pdfoutput' is a single logical")
  
  
  # number of unique antibodies
  numAB <- length(unique(tidydf$Antibody.Name))
  
  # make plot
  if (pdfoutput){
    pdf("RFI_barplots/plotPerAB.pdf")
  }
  
  
  for(i in 1:numAB){
    
    print(
      tidydf %>%
        filter(Antibody.Name == unique(tidydf$Antibody.Name)[i]) %>%    # filter for only sample i
        ggplot(aes_string(y=RFIcol,x="X1")) + 
        geom_bar(stat = "identity") +
        labs(title = unique(tidydf$Antibody.Name)[i], y = "RFI", 
             x = "Sample") +
        theme(plot.title = element_text(hjust = 0.5),
              title = element_text(size=11)) +
        coord_flip()
    )
  }
  
  if (pdfoutput){
    dev.off()
  }
  
} 