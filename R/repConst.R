#' Replicate consistency
#' 
#' Explore replicate consistency using the mean, standard deviation and 
#' coefficient of variation of the RFI values per antibody of each condition 
#' group.
#' 
#' 
#' @param tidydf Tidy dataframe of RFI data.
#' @param reps Vector of replicate condition names as string, to calculate 
#'     consistency of.
#' @param cond_col Name of the column that contains the condition information.
#' @param title Title of the faceted barplot, as string.
#' @param filename Filename of the replicate consistency data (include '.tsv'
#'     at the end).
#' 
#' 
#' @importFrom assertthat assert_that
#' @importFrom magrittr %>%
#' @importFrom rlang !!
#' @import ggplot2
#' 

#' @describeIn plot_repConst Faceted bar plot of the mean +/- standard 
#'     deviation of each condition group, per antibody.
#' @export
plot_repConst <- function(tidydf, reps, cond_col = "Condition", title) {
  
  # check inputs
  assert_that(sum(c(cond_col, "Antibody.Name", "RFI") %in% colnames(tidydf))
              == 3, msg = "Check your 'cond_col' and 'Antibody.Name' & 
                           'RFI' columns exist in 'tidydf'")
  
  assert_that(is.character(reps), msg = "Check that 'reps' is a string vector")
  
  assert_that(is.character(cond_col), length(cond_col) == 1,
              msg = "Check that 'cond_col' is a single string")
  
  assert_that(is.character(title), length(title) == 1,
              msg = "Check that 'title' is a single string")
  
  
  
  tidydf %>%
    dplyr::filter(!!as.name(cond_col) %in% reps) %>%
    dplyr::group_by(!!as.name(cond_col), Antibody.Name) %>%
    dplyr::summarise(mean = mean(RFI), sd = sd(RFI)) %>%
    ggplot(aes(y = mean, x = !!as.name(cond_col))) + 
    geom_bar(stat = "identity", position = "dodge") + 
    geom_errorbar(aes(ymin=mean-sd, ymax = mean+sd, width = 0.5)) + 
    facet_wrap(Antibody.Name~.) + 
    coord_flip() +
    labs(title = title, y = "Mean RFI for each group (error bar=SD)")
  
}

#' @describeIn plot_repConst Write to file the mean, standard deviation and
#'     coefficient of variation of each condition, per antibody.
#' @export 
df_repConst <- function(tidydf, reps, cond_col = "Condition", filename) {
  
  # check inputs
  assert_that(sum(c(cond_col, "Antibody.Name", "RFI") %in% colnames(tidydf))
              == 3, msg = "Check you 'cond_col' and 'Antibody.Name' & 
                           'RFI' columns exist in 'tidydf'")
  
  assert_that(is.character(reps), msg = "Check that 'reps' is a string vector")
  
  assert_that(is.character(filename), length(filename) == 1,
              msg = "Check that 'filename' is a single string")
  
  
  const_df <- tidydf %>%
    dplyr::filter(!!as.name(cond_col) %in% reps) %>%
    dplyr::group_by(!!as.name(cond_col), Antibody.Name) %>%
    dplyr::summarise(Mean = mean(RFI), SD = sd(RFI), CV = SD/Mean) %>%
    
    
  write.table(const_df, file = filename, sep = "\t", quote = FALSE,
              row.names = FALSE)
  
}






