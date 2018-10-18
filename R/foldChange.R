###########################################################
# Functions to calculate fold change between 
# conditions/samples and visualise them for RPPA data
#
# Author: Lucy Liu
###########################################################

# colour palettes
pal <- c("#000000","#004949","#009292","#ff6db6","#ffb6db",
         "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
         "#920000","#924900","#db6d00","#24ff24","#ffff6d")

pal2 <- c("#000000", "#009E73", "#e79f00", "#9ad0f3", "#0072B2", 
          "#D55E00", "#CC79A7", "#F0E442")

##########################################################
# Calculate fold change given data in tidy (long) format
# and df of comparisons
##########################################################

calcFC <- function(
  tidydf,       # df of RFI data in tidy (long) format. There should  
                # NOT be any technical replicates
  comparisons,  # df of all the comparisons to make, using sample 
                # names
  log = FALSE,  # logical indicating whether RFI values should be 
                # log2'ed before calculated foldchange
  ABnames       # optional argument. If given AB names will be merged
){
  
  # check inputs
  assertthat::assert_that(
    sum(c("X1", "RFI", "AB") %in% colnames(tidydf)) == 3,
    msg = "Check 'data' has columns 'X1', 'RFI' and 'AB'")
  
  assertthat::assert_that(is.logical(log), length(log) == 1,
                          msg = "Check 'log' is single logical")
  
  assertthat::assert_that(sum(
    c(comparisons[,1], comparisons[,2]) %in% tidydf$X1) == 2* nrow(comparisons),
    msg = "Check 'comparisons' df uses sample names")
   
  if (! missing(ABnames)){
    assertthat::assert_that(
      sum(c("Antibody.Name","Ab.No.") %in% colnames(ABnames)) == 2,
      dim(ABnames)[2] == 2,
      msg = "Check column names in 'ABnames' dataframe")
  }
  
  # convert to data.frame
  if (sum(class(tidydf) %in% "tbl_df") >= 1){
    tidydf <- as.data.frame(tidydf)
  }
  
  # obtain number of comparisons and ABs
  num_comparisons <- nrow(comparisons)
  num_ABs <- length(unique(tidydf$AB))
  
  # convert to wide format
  wide_df <- tidydf %>%
    select(X1, RFI, AB) %>%
    spread(value = RFI, key = AB)
  
  # convert to matrix
  numeric_mat <- as.matrix(wide_df[,-1])
  rownames(numeric_mat) <- wide_df[,1]
  
  
  # calculate fold changes
  fc_mat <- matrix(nrow = num_comparisons, 
                   ncol = num_ABs)
  
  if (log){
    
    numeric_mat <- log2(numeric_mat)
  }
  
  for (i in 1:num_comparisons){
      
      cond1 <- comparisons[i,1]
      cond2 <- comparisons[i,2]
      
      if (log){
        
        # subtract if log
        fc_mat[i,] <- numeric_mat[rownames(numeric_mat) == cond2,] - 
          numeric_mat[rownames(numeric_mat) == cond1,]
        
      } else {
        
        # divide if raw
        fc_mat[i,] <- numeric_mat[rownames(numeric_mat) == cond2,] / 
          numeric_mat[rownames(numeric_mat) == cond1,]
        
      }
  }
  
  
  # convert back to df
  fc_df <- as.data.frame(fc_mat)
  
  ## add AB names
  colnames(fc_df) <- colnames(wide_df)[-1]
      
  ## add conditions column
  condition <- vector(mode = "character", length = num_comparisons)
      
  for (i in 1:num_comparisons){
    
    condition[i] <- pheno[pheno$Lysate.ID == comparisons[i,2], "Condition"]
    
  }
  
  fc_df$Condition <- condition
  
  # add sample columns
  fc_df$Sample1 <- comparisons[,1]
  fc_df$Sample2 <- comparisons[,2]
      
  # change to long form
  fc_df <- fc_df %>%
    gather(1:num_ABs, key = "AB", value = "FoldChange") 
  
  if (! missing(ABnames)){
    fc_df <- merge(fc_df, ABnames, by.x = "AB", by.y = "Ab.No.",
                all.x = TRUE, all.y = FALSE)
  }
  
  return(fc_df)
      
}

  

####################################
# Plot fold change graphs
####################################
  
  
plot_fc <- function(
  fc_df,               # fold change df
  samples,             # vector of "cond2" (unique/non-control samples)
                       # sample names to show on plot
  logged = FALSE,      # logical indicating whether the data have been logged
  normalised = FALSE,  # logical indicating whether the data have been 
                       # normalised
  flip                 # logical indicating whether the plot should be flipped
){
  
  # check inputs
  assertthat::assert_that(sum(c("Condition","Sample1","Sample2","AB", 
                                "FoldChange") %in% colnames(fc_df)) == 5,
                          msg = "Check 'fc_df' has the correct columns")
  
  assertthat::assert_that(is.character(samples),
                          msg = "Check that 'samples' is of string type")
  
  assertthat::assert_that(is.logical(logged), length(logged) == 1,
                          msg = "Check 'logged' is a single logical")
  
  assertthat::assert_that(is.logical(normalised), length(normalised) == 1,
                          msg = "Check 'normalised' is a single logical")
  
  assertthat::assert_that(is.logical(flip), length(flip) == 1,
                          msg = "Check 'flip' is a single logical")
  
  
  # create titles
  ylab <- "Fold Change"
  
  if (logged){
    
    ylab <- paste("Log2", ylab)
    
  }
  
  if (normalised){
    
    ylab <- paste("Normalised", ylab)
    
  }
  
  # number of comparisions
  num_comparisions <- nrow(comparisons)
  
  # main plot
  gg <- fc_df %>%
    filter(Sample2 %in% samples) %>%
    mutate(SampCond = paste(Sample2, Condition)) %>%
    ggplot(aes(y = FoldChange, x = Antibody.Name, fill = SampCond)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "Fold change per antibody", x = "Antibody", 
         y = ylab) +
    scale_fill_manual(values = pal2[1:num_comparisions],
                      name = "Sample & Condition") +
    theme(plot.title = element_text(hjust = 0.5))
  
  
  if (flip){
    
    gg <- gg +
      coord_flip()
      
      
  } else {
    
    gg <- gg + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }

  return(gg)
  
}
  
  
  
  
  
  
  
  

