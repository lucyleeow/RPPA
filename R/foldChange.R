#############################################
# Functions to calculate fold change between 
# conditions/samples and visualise them for RPPA data
#
# Author: Lucy Liu
#############################################

###########################################################
# Calculate fold change given data in tidy (long) format
# and df of comparisons
###########################################################

calcFC <- function(
  tidydf,       # df of RFI data in tidy (long) format. There should  
                # NOT be any technical replicates
  comparisons,  # df of all the comparisons to make, using sample 
                # names
  log = FALSE   # logical indicating whether RFI values should be 
                # log2'ed before calculated foldchange
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
          numeric_mat[rownames(numeric_mat) == cond1]
        
      } else{
        
        # divide if raw
        fc_mat[i,] <- numeric_mat[rownames(numeric_mat) == cond2,] / 
          numeric_mat[rownames(numeric_mat) == cond1]
        
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
      
  # change to long form
  fc_df <- fc_df %>%
    gather(1:num_ABs, key = "AB", value = "FoldChange") 
  
  
  return(fc_df)
      
  }
  

  
  
  
  
  
  
  
  
  
  
  
  

