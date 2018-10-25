#' Global rank invariant normalisation (GRIN)
#' 
#' Performs global rank invariant normalisation, as described by \href{https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-520}{Pelz et al.}
#' 
#' The code for this function was taken from \href{https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-520}{Pelz et al.} 
#' Some minor changes have been made to allow it to be used for RPPA 
#' data.
#' 
#' @param data Matrix of RFI values with column for each sample.
#' @param width Width (inches) of diagnostic plot.
#' @param height Height (inches) of diagnostic plot.
#' @param pointsize Point size for text of diagnostic plot.
#' @param filetype Diagnostic plot format "png", "wmf", or "postscript"
#' @param ReportName Name for Diagnostic plots.
#' @param count Size of Global Rank-invariant Set to use.
#' @param f Smoother parameter for lowess.


# #############################################
# Main GRSN implementation script.
# #############################################

#' @export
GRSN <- function(data,           # exprSet from affy package or matrix with column for each sample.
                 width=15,       # Width (inches) of diagnostic plot.
                 height=5,       # Height (inches) of diagnostic plot.
                 pointsize=16,   # Point size for text of diagnostic plot.
                 filetype="png", # Diagnostic plot format "png", "wmf", or "postscript".
                 ReportName=NA,  # Name for Diagnostic plots.
                 count=5000,     # Size of Global Rank-invariant Set to use.
                 f=0.25)         # Smoother parameter for lowess.
{
  
  # Changes are marked with LL (Lucy Liu)
  
  # Check the class of the input data.
  if (class(data) == "exprSet")
  {
    rawData <- attr(data, "exprs")
  } else if (class(data) == "matrix")
  {
    rawData <- data
  } else if (class(data)[[1]] == "ExpressionSet")
  {
    rawData <- exprs(data)
  } else
  {
    stop("data parameter is not a valid type!")
  }
  
  ########
  # LL - comment out this test as not relevant to RPPA data (we would hardly ever get
  # RFI over 31) and we know we will give raw, unlogged data.
  ########
  # # Make sure that the input data is not log scale.
  # if (max(rawData) > 31)
  # { # Assume that input data is not log scale.
  isItLogScaled <- FALSE
  # } else
  # { # Assume that input data is log scale and convert it.
  #   isItLogScaled <- TRUE
  #   rawData = 2^rawData
  # }
  
  ############
  ## LL - comment out as affymetrix is not relevant to RPPA
  ############
  
  
  # Find Affymetrix(R) control probe sets by looking for 
  # probe set IDs starting in "AFFY".
  # affyIdx <- grep ("^AFFX", attr(rawData, "dimnames")[[1]])
  
  # Data to normalize.
  # adjust <- max(0, (0.25 - min(rawData)))
  ## LL: change adjust to prior so not to be confused with the other
  ## adjust variable below
  
  prior <- max(0, (0.01 - min(rawData)))
  # If your smallest data point is bigger than 0.0.1, then adjust = 0, 
  # meaning that we do not need to add anything when we log
  
  # If your smallest data point is smaller than 0.01, we would add the difference
  # between the smallest number and 0.01, such that the smallest value to log is 0.01
  
  ### Log your data
  M1 <- log2(rawData + prior)
  
  #####
  # Note M1 is logged data!
  #####
  
  # Get the average of the reference set.
  # Do a trimmed mean to be robust, but eliminate the "artifact" that 
  # shows up when doing median on an odd number of samples.
  
  
  # Trim removes the 0.25% of observations from each end
  # This gives you the trimmed average RFI of each row (AB)
  Mavg <- apply(M1[, ], 1, mean, trim=0.25)
  
  # New method for a global invariant set.
  total <- dim(M1)[[1]]
  # Gives you the total number of rows
  idx <- 1:total
  subSet <- 1:total
  
  #####
  # This part is not relevant
  ####
  
  # Exclude Affy control probe sets from 
  # approximate global rank invaraint set (GRiS).
  # if (length(affyIdx) > 0)
  # {
  #   total <- total - length(affyIdx)
  #   idx <- idx[-affyIdx]
  #   subSet <- subSet[-affyIdx]
  # }
  
  # Calculate number of probe sets to exclude at each iteration.
  discardNumber <- (total - count) / 4
  
  ####################################################
  ### Main iteration loop to get approximate GRiS. ###
  ####################################################
  while (TRUE)
  {
    total <- floor(max(total - discardNumber, count))
    # Get the bigger number between count (number of features you want at the end), 
    # and the total number features - discardnumber. 
    # Note that total updates with each iteration
    M2 <- cbind(apply(M1[idx, ], 2, rank))
    # Apply the function rank to each column
    V2 <- apply(M2, 1, var)
    # Find the variance along each row (feature)
    subSet <- order(V2, decreasing=FALSE)[1:total]  
    # Order according to the variance, from lowest to highest variance
    # Then subset only the desired amount
    idx <- idx[subSet]
    if (total == count) break
    # Get out of the loop when your total ==  count
  }
  invariantIdx <- idx
  
  # Use invariant set to normalize all samples to the average.
  Mnew <- NULL
  x <- Mavg
  # Mavg is the trimmed mean of the raw data across a row
  
  for (b in 1:dim(M1)[[2]])
  {
    if (!is.na(ReportName))
    {
      PrintToFile <- .GraphSetup(width=width, height=height, 
                                pointsize=pointsize, 
                                ReportName=paste(ReportName, b, sep=""),
                                filetype=filetype)
    }
    
    if (PrintToFile)
    {
      # Plot three graphs side-by-side.
      layout(matrix(c(1,2,3), nrow=1, ncol=3, byrow = TRUE))
      par(mar=c(5,5,2,1)) # c(bottom, left, top, right).
    }
    
    ######################################
    # Make M and A values from ALL the data
    ######################################
    
    y <- M1[,b]
    # y is a row of data (feature)
    
    ### M vs. A transformed data.  ###
    M <- y-x
    # x is the trimmed mean of the raw data across a row
    A <- (y+x)/2
    
    #######################################################
    ### Lowess curve based on M vs. A transformed data. ###
    ### Note that you only use the invariant features   ### 
    ### to make your lowess curve                       ###
    #######################################################
    
    curve <- lowess(x=A[invariantIdx], y=M[invariantIdx], f=f)
    
    ### Create evenly space lookup from calibration curve. ###
    aCurve <- curve[[1]]
    mCurve <- curve[[2]]
    steps <- 1000
    sampleMin <- min(A)
    sampleMax <- max(A)
    step <- (sampleMax - sampleMin) / steps
    position <- seq(sampleMin, sampleMax, length=steps + 1)
    adjust <- array(0,c(steps+1))
    count <- length(aCurve)
    
    idxL <- 1
    idxR <- 2
    for (i in 1:(steps + 1))
    {
      while (idxR < count && position[i] > aCurve[idxR])
      {
        idxR <- idxR + 1
      }
      while ((idxL + 1) < idxR && position[i] > aCurve[idxL + 1])
      {
        idxL <- idxL + 1
      }
      while (idxR < count && aCurve[idxL] >= aCurve[idxR])
      {
        idxR <- idxR + 1
      }
      ## these while loops ensure that the location of position[i]
      ## is between aCurve[idxL] and aCurve[idxR]
      if (aCurve[idxL] < aCurve[idxR])
      {
        adjust[i] <- (((mCurve[idxR] - mCurve[idxL])/(aCurve[idxR] - aCurve[idxL]))
                      # This is just slope: rise/run
                      *(position[i] - aCurve[idxL]) + mCurve[idxL])
      }
    }
    ## The above simply calculates the value of M at each position[i]
    
    
    ### Apply lookup to data.  Can be applied to transformed or untransformed data. ###
    yPrime <- y - adjust[(A - sampleMin) / step + 1.5]
    mPrime <- yPrime - x
    
    Mnew <- cbind(Mnew, yPrime)
    if (PrintToFile)
    {
      sampleName <- attr(rawData,"dimnames")[[2]][b]
      
      plot(x=A, y=M, pch=".", col="blue", 
           main= paste("Scatter Plot ", b, " vs. Reference Set", sep=""), 
           sub="", 
           xlab=paste("(log2(", sampleName, ")+log2(ref))/2", sep = ""), 
           ylab=paste("log2(", sampleName, ")-log2(ref)", sep = ""))
      lines(x=c(-5, 20), y=c(0, 0), col="black")
      
      plot(x=A[invariantIdx], y=M[invariantIdx], pch=".", col="blue", 
           main= paste("Invariant Scatter Plot ", b, " vs. Reference Set ",sep=""), 
           sub="", 
           xlab=paste("(log2(", sampleName, ")+log2(ref))/2", sep = ""), 
           ylab=paste("log2(", sampleName, ")-log2(ref)", sep = ""))
      points(x=A[invariantIdx], y=mPrime[invariantIdx], pch=".", col="red")
      lines(x=c(-5, 20), y=c(0, 0), col="black")
      
      lines(x=position, y=adjust, lwd=2, col="green")
      
      plot(x=A, y=mPrime, pch=".", col="red", 
           main= paste("Normalized Scatter Plot ", b, " vs. Reference Set ", sep=""), 
           sub="", 
           xlab=paste("(log2(", sampleName, ")+log2(ref))/2", sep = ""), 
           ylab=paste("log2(", sampleName, ")-log2(ref)", sep = ""))
      lines(x=c(-5, 20), y=c(0, 0), col="black")
    }
    
    if (PrintToFile)
    {
      # Stop outputting to file.
      dev.off()
    }
  }
  
  if (isItLogScaled == FALSE)
  {
    # Convert back to unlogged
    Mnew <- 2^Mnew
  }
  
  # Check the class of the input data.
  if (class(data) == "exprSet")
  {
    # Update expression set values.
    attr(data, "exprs")[,] <- Mnew[,]
  } else if (class(data) == "matrix")
  {
    # Update matrix values, keeping attributes
    data[,] <- Mnew[,]
  } else if (class(data)[[1]] == "ExpressionSet")
  {
    exprs(data) <- Mnew[,]
  }
  
  return(data)
}

#' @keywords internal
GraphSetup <- function(width=6, height=6, pointsize=12, 
                       ReportName="test", filetype="none")
{
  PrintToFile <- FALSE
  if (filetype == "png")
  {
    PrintToFile <- TRUE
    file <- paste(ReportName, ".png", sep="")
    # Write plots as image files.
    png(filename=file, width=96*width, height=96*height, 
        pointsize=pointsize, bg="white")
  } 
  
  if (filetype == "wmf")
  {
    PrintToFile <- TRUE
    file <- paste(ReportName, ".wmf", sep="")
    # Write plots as windows meta files.
    win.metafile(filename=file, width=width, height=height, 
                 pointsize=pointsize)
  } 
  
  if (filetype == "postscript")
  {
    PrintToFile <- TRUE
    file <- paste(ReportName, ".eps", sep="")
    # Write plots as postscript files.
    postscript(file=file, width=width, height=height, 
               pointsize=pointsize, onefile=FALSE, 
               horizontal=FALSE, paper="special", family="Times")
  } 
  
  return(PrintToFile)
}

