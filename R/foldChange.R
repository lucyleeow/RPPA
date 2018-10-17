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
  
  # 
  
  
  num_comparisons <- nrow(comparisons)
  num_Abs <- length(unique(tidydf$AB))
  
  # convert to wide format
  wide_df <- tidydf %>%
    select(X1, RFI, AB) %>%
    spread(value = RFI, key = AB)
  
  # convert to matrix
  numeric_mat <- as.matrix(wide_df[,-1])
  rownames(numeric_mat) <- wide_df[,1]
  
  
  # calculate fold changes
  fc_mat <- matrix(nrow = nrow(comparisons), 
                   ncol = ncol(numeric_mat))
  
  if (log){
    
    numeric_mat <- log2(numeric_mat)
    
    for (i in 1:nrow(comparisons)){
      
      cond1 <- comparisons[i,1]
      cond2 <- comparisons[i,2]
      
      fc_mat[i,] <- numeric_mat[rownames(numeric_mat) == cond2,] - 
        numeric_mat[rownames(numeric_mat) == cond1]
      
    }else{
      
      for (i in 1:nrow(comparisons)){
        
        cond1 <- comparisons[i,1]
        cond2 <- comparisons[i,2]
        
        fc_mat[i,] <- numeric_mat[rownames(numeric_mat) == cond2,] / 
          numeric_mat[rownames(numeric_mat) == cond1]
      
      }
      
      
      
    
    fc_df <- as.data.frame(fc_mat)
      colnames(fc_df) <- colnames(wide_df)[-1]
      
      condition <- vector(mode = "character", length = nrow(comparisons))
      
      for (i in 1:nrow(comparisons)){
        
        condition[i] <- pheno[pheno$Lysate.ID == comparisons[i,2], "Condition"]
        
      }
      
      fc_df$Condition <- condition
      
      fc_df <- fc_df %>%
        gather(1:num_Abs, key = "AB", value = "FoldChange") 
  }
  

  
  
  
  
  
  
  
  
  
  
  
  
  
  
}
