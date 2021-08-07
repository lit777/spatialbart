function_mean <- function(p.t.node, p.tnode, p.i.nodes, rules, ind_sub=NULL){
  
  if(ind_sub==1){
      xpred <- Xpred[-mis.ind,]
      xcut <- lapply(1:dim(xpred)[2], function(t) sort(unique(xpred[,t]))) # e.g. unique values of predictors
    nn <- n.complete
  }else{
    if(ind_sub==0){
      xpred <- Xpred0; xcut <- Xcut0
    }else{
      xpred <- Xpred; xcut <- Xcut
    }
  }
  
  
  loc.pick <- p.t.node[p.tnode,] # which terminal node?
  
  h.lv <- NULL   # location in a given hierarchical level
  h.lv[1] <- loc.pick[2] 
  lv <- NULL # a hierarchical level
  lv[1] <- loc.pick[1]
  ii <- 1
  while(lv[ii]!=1){
    ii <- ii+1
    if(h.lv[ii-1] %% 2 ==0){
      h.lv[ii] <- h.lv[ii-1]/2
      lv[ii] <- loc.pick[1]-ii+1
    }else{
      h.lv[ii] <- (h.lv[ii-1]+1)/2 
      lv[ii] <- loc.pick[1]-ii+1
    }
  }
  h.lv <- rev(h.lv)
  lv <- rev(lv)
  
  temp <- cbind(xpred, 1:dim(xpred)[1])
  if(h.lv[2] %% 2 ==0 ){
    temp.t <- subset(temp, temp[,p.i.nodes[[1]]] >= xcut[[p.i.nodes[[1]]]][rules[[1]]])
    temp.R <- temp.t[,dim(temp.t)[2]]
    temp.L <- setdiff(1:dim(temp)[1],temp.R)
    R <- temp.R
  } else {
    temp.t <- subset(temp, temp[,p.i.nodes[[1]]] < xcut[[p.i.nodes[[1]]]][rules[[1]]])
    temp.L <- temp.t[,dim(temp.t)[2]]
    temp.R <- setdiff(1:dim(temp)[1],temp.L)
    R <- temp.L
  }
  
  if(length(h.lv) > 2){
    for(i in 3:length(h.lv)){
      if(h.lv[i] %% 2 ==0){
        temp <- temp.t
        temp.t <- subset(temp, temp[,p.i.nodes[[lv[i-1]]][h.lv[i-1]]] >= xcut[[p.i.nodes[[lv[i-1]]][h.lv[i-1]]]][rules[[lv[i-1]]][h.lv[i-1]]])  
        temp.R <- temp.t[,dim(temp.t)[2]]
        temp.L <- setdiff(temp[,dim(temp)[2]], temp.R)
        R <- temp.R
      } else {
        temp <- temp.t
        temp.t <- subset(temp, temp[,p.i.nodes[[lv[i-1]]][h.lv[i-1]]] < xcut[[p.i.nodes[[lv[i-1]]][h.lv[i-1]]]][rules[[lv[i-1]]][h.lv[i-1]]])  
        temp.L <- temp.t[,dim(temp.t)[2]]
        temp.R <- setdiff(temp[,dim(temp)[2]], temp.L)
        R <- temp.L
      }
    }
  }
  
  return(list(R=R))
}

# Sample mean parameters
Mean.Parameter <- function(sigma2, sigma_mu, i.nodes, rules,  R, ind=NULL){
  
  if(ind==1){
    xpred <- Xpred[-mis.ind,]
    xcut <- lapply(1:dim(xpred)[2], function(t) sort(unique(xpred[,t]))) # e.g. unique values of predictors
    nn <- n.complete
  }else{
    if(ind==0){
      xpred <- Xpred0; xcut <- Xcut0; nn <- n0
    }else{
      xpred <- Xpred; xcut <- Xcut; nn <- length(R)
    }
  }
  
  t.node <- c(0, 0)
  for(i in 2:length(i.nodes)){
    temp <- grep(999, i.nodes[[i]])
    if(length(temp) > 0){
      t.node <- rbind(t.node, cbind(i, temp))
    }
  }
  
  if(length(t.node)==2){
    Var <- 1/(1/sigma_mu+length(R)/sigma2)
    Mean <- Var * (sum(R)/sigma2)
    mu <- rnorm(1,Mean, sqrt(Var))
    T <- rep(mu, nn)
  }else{
  
  t.node <- t.node[-1,]
  
  
  Resid <- list()  
  for(i in 1:dim(t.node)[1]){
    Resid[[i]] <- function_mean(t.node, i, i.nodes, rules, ind_sub=ind)$R 
  }  
  
  mu <- NULL
  for(i in 1:dim(t.node)[1]){
    Var <- 1/(1/sigma_mu+length(Resid[[i]])/sigma2)
    Mean <- Var * (sum(R[Resid[[i]]])/sigma2)
    mu[i] <- rnorm(1,Mean, sqrt(Var))
  }
  
  T <- rep(0, nn)
  for(i in 1:length(Resid)){
    temp.ind <- Resid[[i]]
    T[temp.ind] <- mu[i]
  }
  }
  
  return(list(T=T, mu=mu))
}
