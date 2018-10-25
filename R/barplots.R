#' Barplots per sample/antibody
#' 
#' Generate a separate barplot of the RFI value for EACH sample or antibody.
#' 
#' @param tidydf Tidy dataframe of RFI values to plot.
#' @param RFIcol Column name of the RFI column as string.
#' @param pdfoutput Single logical indicating whether the plots should be
#'     output to separate pdf file or printed in console.
#' 
#' 
#' @describeIn plotperSample Generate a barplot of the RFI value for each 
#'     antibody, one barplot per sample. Save in a file called 
#'     'plotPerSample.pdf' if pdfoutput=TRUE.
#'     
#' @importFrom assertthat assert_that
#'     
#'  
#' @export
plotperSample <- function(tidydf, RFIcol = "RFI", pdfoutput = TRUE){
  
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

#' @describeIn plotperSample Generate a barplot of the RFI value for each 
#'     sample, one barplot per antibody. Save in a file called 'plotPerAB.pdf'
#'     if pdfoutput = TRUE.
#' @export
plotperAB <- function(tidydf, RFIcol = "RFI", pdfoutput = TRUE){
  
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