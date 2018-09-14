###############
# RPPA Plots
###############

# These functions generate plots for RPPA data


############################
# RFI for each AB - 1 plot
############################

RFIperAB <- function(
  df,            # dataframe to plot
  reps = FALSE   # logical indicated whether there were replicates
) {
  
  
  # check reps input
  
  assertthat::assert_that(
    is.logical(reps), length(reps) == 1, 
    msg = "Incorrect input for reps argument. Should be single logical."
    )
  
  
  # determine column name of RFI
  
  if (reps){
    colRFI <- "RFImean"
  } else {
    colRF <- "RFI"
  }
  
  # make plot
  
    ggplot(df, aes_string(y=colRFI, x="Antibody.Name")) + 
      geom_boxplot() + 
      labs(title = "RFI per AB", y = "RFI", x = "Antibody") +
      coord_flip() +
      theme(plot.title = element_text(hjust = 0.5), title = element_text(size=16),
            axis.text.y = element_text(size = 14))
 
}

#################################
# RFI for each Sample - 1 plot
#################################

RFIperSample <- function(
  df,            # dataframe to plot
  reps = FALSE   # logical indicated whether there were replicates
) {
  
  # check reps input
  
  assertthat::assert_that(
    is.logical(reps), length(reps) == 1, 
    msg = "Incorrect input for reps argument. Should be single logical."
  )
  
  # determine column name of RFI
  
  if (reps){
    colRFI <- "RFImean"
  } else {
    colRF <- "RFI"
  }
  
  # make plot
  
  ggplot(df, aes_string(y=colRFI, x="X1")) + 
    geom_boxplot() + 
    labs(title = "RFI per Sample", y = "RFI", x = "Sample") +
    coord_flip() +
    theme(plot.title = element_text(hjust = 0.5), title = element_text(size=16),
          axis.text.y = element_text(size = 14))
}


###########################################
# RFI for each Sample - 1 plot per sample
###########################################

plotperSample <- function(
  df      # dataframe to plot
  ){
  
  numSample <- length(unique(df$X1))
  
  for(i in 1:numSample){
   
     print(
      df %>%
        filter(X1 == df$X1[i]) %>%    # filter for only sample i
        ggplot(aes(y=RFI,x=Antibody.Name)) + 
        geom_bar(stat = "identity") +
        labs(title = df$X1[i], y = "Average RFI", x = "Antibody") +
        theme(plot.title = element_text(hjust = 0.5),
              title = element_text(size=14)) +
        coord_flip()
    )
  }


}


##################################
# RFI for each AB - 1 plot per AB
##################################

plotperAB <- function(
  df      # dataframe to plot
){
  
  numAB <- length(unique(df$AB))
  
  for(i in 1:numAB){
    
    print(
      df %>%
        filter(Antibody.Name == df$Antibody.Name[i]) %>%    # filter for only sample i
        ggplot(aes(y=RFI,x=X1)) + 
        geom_bar(stat = "identity") +
        labs(title = df$Antibody.Name[i], y = "RFI", 
             x = "Sample") +
        theme(plot.title = element_text(hjust = 0.5),
              title = element_text(size=14)) +
        coord_flip()
    )
  }
  
}










