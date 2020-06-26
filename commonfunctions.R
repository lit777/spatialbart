###################################################
# Check the W matrix
# This code is from CARBayes R package (Duncan Lee)
###################################################

common.Wcheckformat <- function(W){
  #### Check W is a matrix of the correct dimension
  if(!is.matrix(W)) stop("W is not a matrix.", call.=FALSE)
  n <- nrow(W)
  if(ncol(W)!= n) stop("W is not a square matrix.", call.=FALSE)    
  
  #### Check validity of inputed W matrix
  if(sum(is.na(W))>0) stop("W has missing 'NA' values.", call.=FALSE)
  if(!is.numeric(W)) stop("W has non-numeric values.", call.=FALSE)
  if(min(W)<0) stop("W has negative elements.", call.=FALSE)
  if(sum(W!=t(W))>0) stop("W is not symmetric.", call.=FALSE)
  
  melted.W <- melt(t(W))
  melted.W <- melted.W[,c(2,1,3)]
  W.triplet <- subset(melted.W, value > 0)
  names(W.triplet) <- NULL
  W.triplet <- as.matrix(W.triplet)
  n.triplet <- nrow(W.triplet) 
  W.triplet.sum <- tapply(W.triplet[ ,3], W.triplet[ ,1], sum)
  n.neighbours <- tapply(W.triplet[ ,3], W.triplet[ ,1], length)
  
  #### Create the start and finish points for W updating
  W.begfin <- array(NA, c(n, 2))     
  count <- 1
  for(i in 1:n){
    W.begfin[i, ] <- c(count, (count + n.neighbours[i]-1))
    count <- count + n.neighbours[i]
  }
  
  #### Return the critical quantities
  results <- list(W=W, W.triplet=W.triplet, n.triplet=n.triplet, W.triplet.sum=W.triplet.sum, n.neighbours=n.neighbours, W.begfin=W.begfin, n=n)
  return(results)   
}

