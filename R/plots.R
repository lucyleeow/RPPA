#############################################
# RPPA Plots
# These functions generate plots for RPPA data
#
# Author: Lucy Liu
#############################################

############################
# RFI for each AB - 1 plot
############################

RFIperAB <- function(
  df,               # dataframe to plot
  RFIcol = "RFI",   # name of the RFI column
  fillby            # optional argument. Fill by column
) {
  
  if (missing(fillby)){
    
    gg <- ggplot(df, aes_string(y=RFIcol, x="Antibody.Name"))
    
  } else{
    
    gg <- ggplot(df, aes_string(y=RFIcol, x="Antibody.Name", fill=fillby))
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

#################################
# RFI for each Sample - 1 plot
#################################

RFIperSample <- function(
  df,               # dataframe to plot
  RFIcol = "RFI",   # name of the RFI column
  fillby            # optional argument. Fill by column
) {
  
  if (missing(fillby)){
    
    gg <- ggplot(df, aes_string(y=RFIcol, x="X1"))
    
  } else{
    
    gg <- ggplot(df, aes_string(y=RFIcol, x="X1", fill=fillby))
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


###########################################
# RFI for each Sample - 1 plot per sample
###########################################

plotperSample <- function(
  df,               # dataframe to plot
  RFIcol = "RFI"    # name of the RFI column
  ){
  
  # number of unique samples
  numSample <- length(unique(df$X1))
  
  for(i in 1:numSample){
   
     print(
      df %>%
        filter(X1 == unique(df$X1)[i]) %>%    # filter for only sample i
        ggplot(aes_string(y=RFIcol,x="Antibody.Name")) + 
        geom_bar(stat = "identity") +
        labs(title = unique(df$X1)[i], y = "Average RFI", x = "Antibody") +
        theme(plot.title = element_text(hjust = 0.5),
              title = element_text(size=11), 
              axis.text.y = element_text(size = 6)) +
        coord_flip()
    )
  }


}


##################################
# RFI for each AB - 1 plot per AB
##################################

plotperAB <- function(
  df,               # dataframe to plot
  RFIcol = "RFI"    # name of the RFI column
){
  
  # number of unique antibodies
  numAB <- length(unique(df$Antibody.Name))
  
  for(i in 1:numAB){
    
    print(
      df %>%
        filter(Antibody.Name == unique(df$Antibody.Name)[i]) %>%    # filter for only sample i
        ggplot(aes_string(y=RFIcol,x="X1")) + 
        geom_bar(stat = "identity") +
        labs(title = unique(df$Antibody.Name)[i], y = "RFI", 
             x = "Sample") +
        theme(plot.title = element_text(hjust = 0.5),
              title = element_text(size=11)) +
        coord_flip()
    )
  }
  
}










