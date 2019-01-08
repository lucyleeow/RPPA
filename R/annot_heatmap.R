#' Annotated heatmap
#' 
#' Generate an annotated heatmap using pheatmap packaged.
#' 
#' @param mat Numeric matrix where each row is an antibody and each column is
#'     a sample. Column names should be sample names.
#' @param tidydf Tidy dataframe of RFI and annotation data.
#' @param annot1 Name of the column in tidydf containing the annotation data,
#'     as string.
#' @param annot2 Optional argument. Name of the column in tidydf containing 
#'     the second annotation data, as string.
#' @param title Title of the graph as string.
#' @param scale Logical. If TRUE, the heatmap will be scaled across rows 
#'     (antibodies), if FALSE no scaling will be done.
#' 
#' @importFrom assertthat assert_that
#' 
#' @export 
annot_heatmap <- function(mat, tidydf, annot1, annot2, title, scale) {
  
  # check inputs
  assert_that(is.numeric(mat), is.matrix(mat), 
             sum(colnames(mat) %in% tidydf$X1) == ncol(mat),
             msg = "Check 'mat' is a numeric matrix with sample names as
             column names")
  
  assert_that(is.character(annot1), length(annot1) == 1, 
              msg = "Check that 'annot1' is a single string")
  
  if (missing(annot2)) {
    assert_that(sum(c(annot1, "X1") %in% colnames(tidydf)) == 2,
                msg = "Check X1 and 'annot1' columns exist in 'tidydf'")
  } else {
    assert_that(is.character(annot2), length(annot2) == 1, 
                msg = "Check that 'annot2' is a single string")
    
    assert_that(sum(c(annot1, annot2, "X1") %in% colnames(tidydf)) == 3,
                msg = "Check X1, 'annot1' and 'annot2' columns exist in 
                'tidydf'")
  }

  assert_that(is.character(title), length(title) == 1,
              msg = "Check 'title' is a single string")
  
  assert_that(is.logical(scale), length(scale) == 1,
              msg = "Check that 'scale' is a single logical")
  
  
  # make annotation dataframe
  column_annot <- data.frame(tidydf[match(colnames(mat), tidydf$X1), annot1])
  colnames(column_annot) <- annot1
  
  if (! missing(annot2)) {
    
    column_annot[,2] <- tidydf[match(colnames(mat), tidydf$X1), annot2]
    colnames(column_annot)[2] <- annot2
    
  }
  
  rownames(column_annot) <- colnames(mat)
  
  
  # make colour vector
  
  ## annot1
  num_annot1 <- length(unique(column_annot[,1]))
  
  assert_that(num_annot1 >=2, msg = "Check that there are 2 or more conditions
              for 'annot1'")
  
  if (num_annot1 == 2) {
    
    annot1_colours <- c("#1B9E77", "#D95F02" )
    names(annot1_colours) <- unique(column_annot[,1])
    
  } else {
    
    annot1_colours <- RColorBrewer::brewer.pal(num_annot1, "Dark2")
    names(annot1_colours) <- unique(column_annot[,1])
    
  }
  
  
  ## annot2
  if (! missing(annot2)) {
    
    num_annot2 <- length(unique(column_annot[,2]))
    
    assert_that(num_annot2 >=2, msg = "Check that there are 2 or more conditions
              for 'annot1'")
    
    if (num_annot2 == 2) {
      
      annot2_colours <- c("#A6CEE3", "#1F78B4" )
      names(annot2_colours) <- unique(column_annot[,2])
      
    } else {
      annot2_colours <- RColorBrewer::brewer.pal(num_annot2, "Paired")
      names(annot2_colours) <- unique(column_annot[,2])
    }
    
  }
  
  annot_colours <- list(annot1_colours)
  names(annot_colours) <- annot1
  
  if (! missing(annot2)) {
    
    annot_colours[[2]] <- annot2_colours
    names(annot_colours)[2] <- annot2
  }
  
  
  if (scale) {
    scale <- "row"
  } else {
    
    scale <- "none"
  }
  
  
  pheatmap::pheatmap(mat[!rownames(mat) == "Ab-41",], fontsize = 8, 
           annotation_col = column_annot, annotation_colors = annot_colours,
           main = title, scale = scale)
  
  
}