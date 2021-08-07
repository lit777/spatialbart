#################
# Simulation Data
#################

# Load reauired libraries
library(MASS)
library(mnormt)
library(MCMCpack)
library(rootSolve)
library(R.matlab)
library(truncnorm)
library(data.table)
library(Rcpp)

  # Square lattice region
  x.easting <- 1:8
  x.northing <- 1:8
  Grid <- expand.grid(x.easting, x.northing)
  K <- nrow(Grid)
  
  # Define spatial connections (see Kim (2020+) for more details)
  distance <- as.matrix(dist(Grid))
  Distance <- distance
  diag(Distance) <- 1
  W <-array(0, c(K,K))
  W[distance==1] <-0.1
  W.hole <- read.csv("W_sim2.csv")
  W.hole <- W.hole[-c(1),-c(1:2)]
  W.hole <- ifelse(is.na(W.hole), 0.1, 1000)
  colnames(W.hole) <- NULL
  rownames(W.hole) <- NULL
  diag(W.hole) <- 0
  wind_mat <- W.hole
  wind_mat <- ifelse(wind_mat==1000, 0, 1)
  
  
  #### Generate the covariates and response data
  divide.ind <- c(rep(1,4),rep(0,4),rep(1,4),rep(0,4),rep(1,4),rep(0,4),rep(1,5),rep(0,3),rep(1,5),rep(0,3),rep(1,5),rep(0,3),rep(1,8),rep(1,8))
  x1 <- divide.ind*rnorm(K, 2.5, 0.25) + (1-divide.ind)*rnorm(K, 1, 0.25)
  x2 <- rnorm(K, 1-0.1*(Grid[,2]-10)^2, 0.1)
  x.com <- rmnorm(K, c(-2, 2), matrix(c(0.3,0.1,0.1,0.3), 2, 2))
  x3 <- x.com[,1]
  x4 <- x.com[,2]
  x5 <- rnorm(K, -1, 0.5) 
  
  theta <- rnorm(K, mean=rep(0,K), sd=0.1)
  diag(Distance) <- 0
  phi <- divide.ind*rmnorm(1, rep(2,K), 0.5*exp(-0.1 * Distance))+(1-divide.ind)*rmnorm(1, rep(-2,K), 0.5*exp(-0.1 * Distance))
  phi <- phi-mean(phi)
  Y <- 1+1*x1 + 1.5*x2 + 1.5*x3 - 1*x4 - 1*x5 + theta + phi
  Y.true <- Y
  
  # randomly selected missing locations
  missing <- c(1:64)[-c(1,3,5,8,9,12,14,16,18,21,27,30,32,35,39,42,44,46,48,49,52,55,57,59,64)]
  Y[missing] <- NA

  P <- 40      # Num. of predictors
  n <- 64
  cov <- list()
  for(i in 1:35){
    cov[[i]] <- rnorm(n,0,1) # draw x independently from a normal distribution
  }
  Xpred1 <- do.call(cbind, cov)
  Xpred <- cbind(x1,x2,x3,x4,x5, Xpred1)
  
  ####################
  # Set the BART model
  ####################
  
  # Import files for the spatial random effect
  sourceCpp("CARBayes.cpp") # this C++ code from CARBayes
  source("commonfunctions.R")
  
  # Import files for the BART model
  source("GROW.R")
  source("PRUNE.R")
  source("CHANGE.R")
  source("Common.R")
  source("Prediction.R")
  
  P <- dim(Xpred)[2]  # Num. of predictors
  Xcut <- lapply(1:dim(Xpred)[2], function(t) sort(unique(Xpred[,t]))) # unique values of the predictors
  shift <- mean(Y, na.rm=TRUE)
  Y <- Y - shift      # Shifting the mean of Y
  
  A.WEIGHT <- NULL    # uncertainty parameter 
  n.full <- length(Y)  # Num. of the locations
  n.complete <- length(which(!is.na(Y)))  # Num. of the locations with observation
  m <- 100 # Num. of the trees in the BART model (default:100)
  n.iter <- 15000 # Num. of MCMC iterations
  
  # Create empty objects
  Tree1 <- matrix(0,nrow=n.complete, ncol=m)  # tree for complete model
  Effect <- matrix(nrow=n.full, ncol=n.iter)
  TREE1 <- list()            # Predicted mean Values
  TREE1[[1]] <- list(T=matrix(nrow=n.complete,ncol=m))
  MU1 <- list()              # Mean parameters associated with terminal nodes
  MU1[[1]] <- list()
  I.nodes1 <- list()         # Place holder of internal nodes (predcitors) + terminal nodes (0's)
  I.nodes1[[1]] <- list()
  Rules1 <- list()           # Place holder of internal nodes (rules)
  Rules1[[1]] <- list()
  Prob <- list()
  pr <- matrix(nrow=P, ncol=1)
  Prob[[1]] <- pr 
  Sigma2 <- NULL            # Variance parameter
  ind <- matrix(nrow=1, ncol=P)
    
  ### New place holder for the next iteration
  TREE1[[2]] <- TREE1[[1]]
  TREE1[[2]] <- TREE1[[1]]
  MU1[[2]] <- MU1[[1]]
  I.nodes1[[2]] <- I.nodes1[[1]]
  Rules1[[2]] <- Rules1[[1]]
  Prob[[2]] <-list(pr=Prob[[1]])
  
  
  # Set (hyper-)parameters and initial values
  prior.tau2 <- c(1, 0.01)
  tau2.posterior.shape <- prior.tau2[1] + 0.5*n.full
  proposal.sd.rho <- 0.02
  rho <- 0.5
  spatial <- rnorm(length(Y), 0, 1)
  tau2 <- 0.1
  R <- Y  # Initial values of R
  
  # Missing index & Obj for imputed outcome
  mis.ind <- which(is.na(Y))
  Y.da <- matrix(nrow=length(mis.ind), ncol=n.iter)
  Y.da[,1] <- rnorm(length(mis.ind), mean(Y, na.rm=TRUE), 1)
  Y[mis.ind] <- Y.da[,1]
  count <- 0
  sigma2 <- 1
  Sigma2[1] <- sigma2
  sigma_mu <- max((min(Y)/(-2*sqrt(m)))^2, (max(Y)/(2*sqrt(m)))^2) # sigma_mu based on min/max of Y
  
  # Set 'nu' and 'q' parameters (Chipman et al. 2010)
  nu <- 3                   # default setting (nu, q) = (3, 0.90) from Chipman et al. 2010
  f <- function(lambda) invgamma::qinvgamma(0.90, nu/2, rate = lambda*nu/2, lower.tail = TRUE, log.p = FALSE) - sqrt(sigma2)
  lambda <- uniroot.all(f, c(0.1^5,10))
  
  dir.alpha <- 1 # alpha parameter in the prior for selection probabilities
  post.dir.alpha <- rep(dir.alpha, P)
  alpha <- 0.95             # alpha (1+depth)^{-beta} where depth=0,1,2,...
  beta <- 2                 # default setting (alpha, beta) = (0.95, 2)
  prop.prob <- rdirichlet(1, post.dir.alpha)
  Prob[[2]]$pr[,1] <- prop.prob
  
  # Prob. of selecting an alteration (GROW, PRUNE and CHANGE)
  p.grow <- 0.28            # Prob. of GROW
  p.prune <- 0.28           # Prob. of PRUNE
  p.change <- 0.44          # Prob. of CHANGE
  
  ### Initial Setup for the "m" trees (depth = 0 and 1)
  for(t in 1:m){
    i.nodes1 <- list()
    for(i in 1:2){
      i.nodes1[[i]] <- rep(NA, 2^(i-1))
    }
    i.nodes1[[1]] <- 999
    I.nodes1[[1]][[t]] <- i.nodes1
    
    rules1 <- list()
    for(i in 1:2){
      rules1[[i]] <- rep(NA, 2^(i-1))
    }
    Rules1[[1]][[t]] <- rules1
  }
  
  # Weight matrix
  a.weight <- 0.95 # initial value for a_w parameter in the weight
  
  diag(distance) <- NA 
  W2 <- 1/distance
  diag(W2) <- 0
  W2 <- ifelse(W2 == Inf, 0, W2)
  
  W2.wind <- wind_mat*a.weight*W2+(1-wind_mat)*(1-a.weight)*W2
  W2.wind <- W2.wind[-mis.ind, -mis.ind]
  rownames(W2.wind) <- 1:(n.full-length(mis.ind))
  colnames(W2.wind) <- 1:(n.full-length(mis.ind))
  W2.post <- common.Wcheckformat(W2.wind)
  
  W2.wind.full <- wind_mat*a.weight*W2+(1-wind_mat)*(1-a.weight)*W2
  W2.post.full <- common.Wcheckformat(W2.wind.full)
  
  det.Q <- 0
  Wstar <- 0
  Wstar.eigen <- 0
  Wstar.val <- 0
  
  #####################################################################
  ############# Run main MCMC #########################################
  #####################################################################
  
  for(j in 2:n.iter){
    count <- count + 1
    
    for(t in 1:m){
      
      R <- Y[-mis.ind] - rowSums(Tree1[,-t]) - spatial[-mis.ind]
      i.nodes <- I.nodes1[[1]][[t]]
      rules <- Rules1[[1]][[t]]
      
      ### Find the depth of the tree (0 or 1 or 2)
      tree.length <- length(which(!is.na(unlist(rules))))
      if(tree.length == 0){  # tree has no node yet
        grow.step <- GROW.first(sigma2=Sigma2[j-1], sigma_mu=sigma_mu, i.nodes=I.nodes1[[1]][[t]],rules=Rules1[[1]][[t]], R=R, prop.prob=prop.prob, ind=1)
        I.nodes1[[2]][[t]] <- grow.step$i.nodes
        Rules1[[2]][[t]] <- grow.step$rules
      } else {
        step <- sample(1:3, 1, prob=c(p.grow, p.prune, p.change))  # Pick a step
        if(step==1){  # GROW step
          grow.step <- GROW(sigma2=Sigma2[j-1], sigma_mu=sigma_mu, i.nodes=I.nodes1[[1]][[t]],rules=Rules1[[1]][[t]], R=R, prop.prob=prop.prob, ind=1)
          I.nodes1[[2]][[t]] <- grow.step$i.nodes
          Rules1[[2]][[t]] <- grow.step$rules
        }
        if(step==2){   # PRUNE step
          prune.step <- PRUNE(sigma2=Sigma2[j-1], sigma_mu=sigma_mu, i.nodes=I.nodes1[[1]][[t]],rules=Rules1[[1]][[t]], R=R, prop.prob=prop.prob, ind=1)
          I.nodes1[[2]][[t]] <- prune.step$i.nodes
          Rules1[[2]][[t]] <- prune.step$rules
        }
        if(step==3){   # CHANGE step
          change.step <- CHANGE(sigma2=Sigma2[j-1], sigma_mu=sigma_mu, i.nodes=I.nodes1[[1]][[t]],rules=Rules1[[1]][[t]], R=R, prop.prob=prop.prob, ind=1)
          I.nodes1[[2]][[t]] <- change.step$i.nodes
          Rules1[[2]][[t]] <- change.step$rules
        }
      }
      
      Mean <- Mean.Parameter(Sigma2[j-1], sigma_mu, i.nodes=I.nodes1[[2]][[t]], rules=Rules1[[2]][[t]],  R, ind=1)
      T <- Mean$T
      mu <- Mean$mu
      
      eval(parse(text=(paste("TREE1[[",2,"]]$T[,",t,"]=T",sep=""))))
      eval(parse(text=(paste("Tree1[,",t,"]=T",sep=""))))
      eval(parse(text=(paste("MU1[[",2,"]]$mu",t,"=mu",sep=""))))
    }
    
    # Sample variance parameter 
    Sigma2[j] <- 1 / rgamma(1, 1+0.5*n.complete, scale=1/(0.01+0.5*sum((Y[-mis.ind]-rowSums(Tree1)-spatial[-mis.ind])^2)))
    
    #######################################
    # Start update the spatial random effect
    ####################################### 
    offset <- (Y[-mis.ind]-rowSums(Tree1))
    spatial[-mis.ind] <- gaussiancarupdate(Wtriplet=W2.post$W.triplet, Wbegfin=W2.post$W.begfin, W2.post$W.triplet.sum, nsites=length(Y[-mis.ind]), phi=spatial[-mis.ind], tau2=tau2, rho=rho, nu2=Sigma2[j], offset=offset)
    
    temp2 <- quadform(as.matrix(W2.post$W.triplet), W2.post$W.triplet.sum, W2.post$n.triplet, length(Y[-mis.ind]), spatial[-mis.ind], spatial[-mis.ind], rho)
    tau2.posterior.scale <- temp2 + prior.tau2[2] 
    tau2 <- 1 / rgamma(1, tau2.posterior.shape, scale=(1/tau2.posterior.scale))
    
    Wstar <- diag(apply(W2.wind,1,sum)) - W2.wind
    Wstar.eigen <- eigen(Wstar)
    Wstar.val <- Wstar.eigen$values
    det.Q <- 0.5 * sum(log((rho * Wstar.val + (1-rho))))    
    
    # update rho parameter
    proposal.rho <- rtruncnorm(n=1, a=0, b=1, mean=rho, sd=proposal.sd.rho)  
    temp3 <- quadform(as.matrix(W2.post$W.triplet), W2.post$W.triplet.sum, W2.post$n.triplet, length(Y[-mis.ind]), spatial[-mis.ind], spatial[-mis.ind], proposal.rho)
    tau3.posterior.scale <- temp3 + prior.tau2[2]
    tau3 <- 1 / rgamma(1, tau2.posterior.shape, scale=(1/tau3.posterior.scale))
    det.Q.proposal <- 0.5 * sum(log((proposal.rho * Wstar.val + (1-proposal.rho))))              
    logprob.current <- det.Q - temp2 / tau2
    logprob.proposal <- det.Q.proposal - temp3 / tau3
    hastings <- log(dtruncnorm(x=rho, a=0, b=1, mean=proposal.rho, sd=proposal.sd.rho)) - log(dtruncnorm(x=proposal.rho, a=0, b=1, mean=rho, sd=proposal.sd.rho)) 
    prob <- exp(logprob.proposal - logprob.current + hastings)
    
    #### Accept or reject the proposal
    if(prob > runif(1))
    {
      rho <- proposal.rho
      det.Q <- det.Q.proposal
    }              
    
    
    for(i in seq_along(mis.ind)){
      spatial[mis.ind[i]] <- rnorm(1, rho*sum(W2.wind.full[mis.ind[i], -mis.ind[i]]*spatial[-mis.ind[i]])/(rho*W2.post.full$W.triplet.sum[mis.ind[i]]+1-rho),  sqrt(tau2 /(rho*W2.post.full$W.triplet.sum[mis.ind[i]]+1-rho) ) )
    }

    # update a_w parameter in the weights  
        proposal.a.weight <- min(rbeta(1, a.weight/(1-a.weight)*2, 2), 1-0.1^7)
    
        W2.wind.proposal <- wind_mat*proposal.a.weight*W2+(1-wind_mat)*(1-proposal.a.weight)*W2
        W2.wind.proposal <- W2.wind.proposal[-mis.ind, -mis.ind]
        rownames(W2.wind.proposal) <- 1:(n.full-length(mis.ind))
        colnames(W2.wind.proposal) <- 1:(n.full-length(mis.ind))
        W2.post.proposal <- common.Wcheckformat(W2.wind.proposal)
    
        temp2 <- quadform(as.matrix(W2.post$W.triplet), W2.post$W.triplet.sum, W2.post$n.triplet, length(Y[-mis.ind]), spatial[-mis.ind], spatial[-mis.ind], rho)
        temp3 <- quadform(as.matrix(W2.post.proposal$W.triplet), W2.post.proposal$W.triplet.sum, W2.post.proposal$n.triplet, length(Y[-mis.ind]), spatial[-mis.ind], spatial[-mis.ind], rho)
        tau3.posterior.scale <- temp3 + prior.tau2[2]
        tau3 <- 1 / rgamma(1, tau2.posterior.shape, scale=(1/tau3.posterior.scale))
        
        Wstar.proposal <- diag(apply(W2.wind.proposal,1,sum)) - W2.wind.proposal
        Wstar.eigen.proposal <- eigen(Wstar.proposal)
        Wstar.val.proposal <- Wstar.eigen.proposal$values
  
        det.Q.proposal <- 0.5 * sum(log((rho * Wstar.val.proposal + (1-rho))))              
        logprob.current <- det.Q - temp2 / tau2
        logprob.proposal <- det.Q.proposal - temp3 / tau3
        hastings <- log(dbeta(proposal.a.weight, 2,1))-log(dbeta(a.weight, 2, 1))+log(dbeta(a.weight, proposal.a.weight/(1-proposal.a.weight)*2, 2)) - log(dbeta(x=proposal.a.weight, a.weight/(1-a.weight)*2, 2)) 
        prob <- exp(logprob.proposal - logprob.current + hastings)
    
    #### Accept or reject the proposal
        if(prob > runif(1))
        {
          a.weight <- proposal.a.weight
          det.Q <- det.Q.proposal
          W2.wind <- W2.wind.proposal
          W2.post <- W2.post.proposal
          Wstar <- Wstar.proposal
          Wstar.eigen <- Wstar.eigen.proposal
          Wstar.val <- Wstar.val.proposal
        }              
    A.WEIGHT[j] <- a.weight
       
    # Update alpha parameter in the Dirichlet prior for the selection probabilities
    add <- as.numeric(table(factor(unlist(I.nodes1[[2]]), levels = 1:P)))
    if(j<n.iter/10){
      post.dir.alpha <- rep(1, P) + add 
    }else{
      p.dir.alpha <- max(rnorm(1, dir.alpha, 0.1), 0.1^10)
      SumS <-  log(ifelse(Prob[[2]]$pr<0.1^300, 0.1^300, Prob[[2]]$pr))
      dir_lik.p <- sum(SumS*(rep(p.dir.alpha/P, P)-1)) + lgamma(sum(rep(p.dir.alpha/P, P)))-sum(lgamma(rep(p.dir.alpha/P, P)))
      dir_lik <- sum(SumS*(rep(dir.alpha/P, P)-1)) + lgamma(sum(rep(dir.alpha/P, P)))-sum(lgamma(rep(dir.alpha/P, P)))
      ratio <- dir_lik.p + log((p.dir.alpha/(p.dir.alpha+P))^(0.5-1)*(P/(p.dir.alpha+P))^(1-1)*abs(1/(p.dir.alpha+P)-p.dir.alpha/(p.dir.alpha+P)^2 ) ) + dnorm(dir.alpha, p.dir.alpha,0.1, log=TRUE) - dir_lik - log((dir.alpha/(dir.alpha+P))^(0.5-1)*(P/(dir.alpha+P))^(1-1)*abs(1/(dir.alpha+P)-dir.alpha/(dir.alpha+P)^2 ) ) - dnorm(p.dir.alpha, dir.alpha,0.1, log=TRUE)
      if (log(runif(1))<ratio) {
        dir.alpha <- p.dir.alpha
      }
      post.dir.alpha <- rep(dir.alpha/P, P) + add 
    }
    
    # update selection prob.
    prop.prob <- rdirichlet(1, post.dir.alpha)
    Prob[[2]]$pr[,1] <- prop.prob
    
    # Store the updates
    TREE1[[1]] <- TREE1[[2]]
    MU1[[1]] <- MU1[[2]]
    I.nodes1[[1]] <- I.nodes1[[2]]
    Rules1[[1]] <- Rules1[[2]]
    Prob[[1]] <-Prob[[2]]
    
    # Predict the outcomes based on the current parameters    
    Tree11 <- matrix(0,nrow=n.full, ncol=m)
    for(t in 1:m){
      i.nodes <- I.nodes1[[2]][[t]]
      rules <- Rules1[[2]][[t]]
      eval(parse(text=(paste("mu <- MU1[[",2,"]]$mu",t,sep=""))))
      
      T <- Mean.Parameter_pred(sigma2=Sigma2[j], sigma_mu=sigma_mu, i.nodes=i.nodes, rules=rules,  mu=mu, ind=3, xpred=Xpred)$T
      eval(parse(text=(paste("Tree11[,",t,"]=T",sep=""))))
    }
    
    # Data augmentation
    Y.da[,count] <- rnorm(length(mis.ind), c(rowSums(Tree11)+spatial)[mis.ind], sqrt(Sigma2[j]))
    Y[mis.ind] <- Y.da[,count]
    ind.temp<-ifelse(add > 0, 1, 0)
    ind <- rbind(ind, ind.temp)
  }
  
  # Post-processing; burn-in the first 7000 iterations and shift the estimates by 'shift'.
  mat <- matrix(nrow=n.full, ncol=5000)
  for(i in 1:5000){
    Y[mis.ind] <- Y.da[,i+9999]
    mat[,i] <- Y+shift
  }
  ind <- ind[10001:15000,]
  mat <- rowMeans(mat)
  
