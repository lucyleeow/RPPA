#' Time series graph
#' 
#' Creates faceted line graphs over time points for each antibody. Each graph
#' can also be coloured by a condition. One line graph will be produced per 
#' antibody and 12 graphs will be displayed per page. The graphs will be 
#' generated in the current graphics device.
#' 
#' @param tidydf Tidy dataframe of RFI values, containing ONLY the 
#'     conditions to plot.
#' @param xcol Column name of the desired x-axis variable, as string.
#' @param xlab x-axis label, as string.
#' @param cond_col Optional argument. Column name of the desired colour 
#'     aesthetic, as string. If given will colour the lines by the condition.
#' @param log Logical indicating whether the RFI data should be logged 
#'     (base 2). A small number (0.00001) is added to prevent taking the log
#'     of zero.
#'     
#' @importFrom assertthat assert_that
#' @importFrom magrittr %>%
#' @import ggplot2
#' 
#' 
#' @export
time_plot <- function(tidydf, xcol, xlab, cond_col, log) {
  
  # check inputs
  assert_that(sum(c("RFI", xcol, "Antibody.Name") %in% colnames(tidydf)) == 3,
              msg = "Check 'RFI', 'Antibody.Name' and your 
              'xcol' columns exist in 'tidydf'")
  
  assert_that(is.character(xcol), length(xcol) == 1,
              msg = "Check 'xcol' is a single string")
  
  assert_that(is.character(xlab), length(xlab) == 1,
              msg = "Check 'xlab' is a single string")
  
  assert_that(is.logical(log), length(log) == 1,
              msg = "Check 'log' is a single logical")
  
  if (! missing(cond_col)) {
    assert_that(is.character(cond_col), length(cond_col) == 1,
                sum(colnames(tidydf) == cond_col) == 1,
                msg = "Check 'cond_col' is a single string and exists in
                'tidydf'")
  }
  
  # calculate number of pages required
  unique_ABs <- unique(tidydf$Antibody.Name)
  num_pages <- ceiling(length(unique_ABs)/12)
  
  # log RFI
  if(log) {
    
    tidydf[["RFI"]] <- log2(tidydf[["RFI"]] + 0.00001)
    ylab = "log2 RFI"
    
  } else {
    
    ylab = "RFI"
    
  }
  
  # make plot
  if (! missing(cond_col)) {
    
    num_conds <- length(unique(tidydf[[cond_col]]))
    
    # RColorBrewer is annoying and gives error if n<2
    if (length(num_conds) <= 2) {
      pal <- c("#1B9E77", "#D95F02" )
    } else {
      pal <- RColorBrewer::brewer.pal(num_conds, "Dark2")
    }
      
    gg <- ggplot(tidydf, aes_string(y = "RFI", x = xcol, 
                                      colour = cond_col)) + 
      geom_point() + 
      stat_summary(aes_string(group = cond_col), fun.y = mean, geom = 'line') + 
      scale_colour_manual(name = "Condition", values = pal)
      
    } else {
      
      gg <- ggplot(tidydf, aes_string(y = "RFI", x = xcol)) + 
        geom_point() + 
        stat_summary(aes_string(group = 1), fun.y = mean, geom = 'line') + 
        scale_colour_manual(name = "Condition")
      
    }
  
  for (i in 1:num_pages) {
    
    print(
      
      gg +
        labs(y = ylab, x = xlab) +
        ggforce::facet_wrap_paginate( ~ Antibody.Name, nrow = 4, ncol = 3, 
                                      page = i) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 40, hjust = 1))
      
    )
    
    # to make sure each plot has it's own page
    cat('\r\n')
    
  }
  
  
}