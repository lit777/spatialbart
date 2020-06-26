# Function of first grow
GROW.first <- function(sigma2, sigma_mu,i.nodes, rules,  R,  prop.prob){

        xpred <- Xpred; xcut <- Xcut
  
        prop.pred <- sample(1:P, 1, replace=FALSE,  prob = prop.prob) # pick a predictor
        prop.rule <- sample(2:length(xcut[[prop.pred]]), 1)
        
        temp <- cbind(xpred, 1:dim(xpred)[1])
        temp.t <- subset(temp, temp[,prop.pred] < xcut[[prop.pred]][prop.rule])
        R.L <- temp.t[,dim(temp.t)[2]]           # Observations < prop.rule
        R.R <- setdiff(1:dim(temp)[1],R.L)    # Observations >= prop.rule
        
        # New tree structure
        prop.i.nodes <- i.nodes
        prop.i.nodes[[1]] <- prop.pred
        prop.i.nodes[[2]] <- c(999,999)
        
        prop.rules <- rules
        prop.rules[[1]] <- prop.rule
        
        # Transition ratio (log scale)
        TRANS <- log(p.prune) + log(1) - log(max(prop.prob[prop.pred],0)) + log(length(xcut[[prop.pred]])-1) - log(p.grow) - log(1)
        
        # Likelihood ratio (log scale)
        nlL <- length(R.L)
        nlR <- length(R.R)
    
        LH <- log(sqrt(sigma2*(sigma2+(nlL+nlR)*sigma_mu))/sqrt((sigma2+nlL*sigma_mu)*(sigma2+nlR*sigma_mu))) + (sigma_mu / (2*sigma2) * ( (sum(R[R.L]))^2/(sigma2+nlL*sigma_mu) +(sum(R[R.R]))^2/(sigma2+nlR*sigma_mu) - (sum(R))^2/(sigma2+(nlR+nlL)*sigma_mu) ) )
    
        # Structure ratio (log scale)
        d <- 0
        STR <- log(alpha) + 2*log((1- alpha / (2 + d)^beta )) - log((1 + d)^beta - alpha) + log(max(prop.prob[prop.pred],0)) - log(length(xcut[[prop.pred]])-1)
    
        r <- TRANS+LH+STR
        if(r > log(runif(1))){
            i.nodes <- prop.i.nodes;
            i.nodes[[length(i.nodes)+1]] <- rep(NA, 2^(length(i.nodes)))
            rules <- prop.rules;
            rules[[length(rules)+1]] <- rep(NA, 2^(length(rules)))
        }
        return(list(i.nodes=i.nodes, rules=rules))
}


# Function of grow alteration
GROW <- function(sigma2, sigma_mu, i.nodes, rules,  R,  prop.prob){

    xpred <- Xpred; xcut <- Xcut
  
    t.node <- c(0, 0)
    for(i in 2:length(i.nodes)){
        temp <- grep(999, i.nodes[[i]])
        if(length(temp) > 0){
            t.node <- rbind(t.node, cbind(i, temp))
        }
    }
    t.node <- t.node[-1,]  # Level and Location of each terminal node
  
    num.t.node <- dim(t.node)[1]
    prop.tnode <- sample(1:num.t.node, 1, replace=FALSE, prob=rep(1/num.t.node, num.t.node)) # pick a terminal node
    
    lv <- NULL # hierarchical level of the selected terminal node
    lv[1] <- t.node[prop.tnode,1]
    h.lv <- NULL   # location in a given hierarchical level
    h.lv[1] <- t.node[prop.tnode,2]
    ii <- 1
    while(lv[ii]!=1){
        ii <- ii+1
        if(h.lv[ii-1] %% 2 ==0){
            h.lv[ii] <- h.lv[ii-1]/2
            lv[ii] <- t.node[prop.tnode,1]-ii+1
        }else{
            h.lv[ii] <- (h.lv[ii-1]+1)/2
            lv[ii] <- t.node[prop.tnode,1]-ii+1
        }
    }
    h.lv <- rev(h.lv)
    lv <- rev(lv)

    temp <- cbind(xpred, 1:dim(xpred)[1])
    if(h.lv[2] %% 2 ==0 ){
        temp.t <- subset(temp, temp[,i.nodes[[1]]] >= xcut[[i.nodes[[1]]]][rules[[1]]])
        temp.R <- temp.t[,dim(temp.t)[2]]
        temp.L <- setdiff(1:dim(temp)[1],temp.R)
    } else {
        temp.t <- subset(temp, temp[,i.nodes[[1]]] < xcut[[i.nodes[[1]]]][rules[[1]]])
        temp.L <- temp.t[,dim(temp.t)[2]]
        temp.R <- setdiff(1:dim(temp)[1],temp.L)
    }
    
    if(length(h.lv) > 2){
        for(i in 3:length(h.lv)){
            if(h.lv[i] %% 2 ==0){
                temp <- temp.t
                temp.t <- subset(temp, temp[,i.nodes[[lv[i-1]]][h.lv[i-1]]] >= xcut[[i.nodes[[lv[i-1]]][h.lv[i-1]]]][rules[[lv[i-1]]][h.lv[i-1]]])
                temp.R <- temp.t[,dim(temp.t)[2]]
                temp.L <- setdiff(temp[,dim(temp)[2]], temp.R)
            } else {
                temp <- temp.t
                temp.t <- subset(temp, temp[,i.nodes[[lv[i-1]]][h.lv[i-1]]] < xcut[[i.nodes[[lv[i-1]]][h.lv[i-1]]]][rules[[lv[i-1]]][h.lv[i-1]]])
                temp.L <- temp.t[,dim(temp.t)[2]]
                temp.R <- setdiff(temp[,dim(temp)[2]], temp.L)
            }
        }
    }

    enough.unique <- which(sapply(1:P, function(x) length(unique(temp.t[,x]))) >= 2) # covariates having enough unique values (2 or more)
    if(length(enough.unique)==0){
        return(list(i.nodes=i.nodes, rules=rules, prob=prop.prob)) # return the current tree if there is no covariate with enough unique values
    }
    prop.pred <- sample((1:P)[enough.unique], 1, replace=FALSE,  prob = prop.prob[enough.unique]) # pick a predictor
    
    unique.len <- length(unique(temp.t[,prop.pred])) # Num. of unique values
    prop.rule <- sort(unique(temp.t[,prop.pred]))[sample(1:(unique.len-1), 1, prob = rep(1/(unique.len-1), (unique.len-1)))+1]
    prop.rule <- which(xcut[[prop.pred]]==prop.rule)

    temp <- temp.t
    temp.t <- subset(temp, temp[,prop.pred] < xcut[[prop.pred]][prop.rule])
    R.L <- temp.t[,dim(temp.t)[2]]
    R.R <- setdiff(temp[,dim(temp)[2]], R.L)
    
    # new tree structure
    prop.i.nodes <- i.nodes
    prop.i.nodes[[t.node[prop.tnode,1]]][t.node[prop.tnode,2]] <- prop.pred
    prop.i.nodes[[t.node[prop.tnode,1]+1]][t.node[prop.tnode,2]*2-1] <- 999
    prop.i.nodes[[t.node[prop.tnode,1]+1]][t.node[prop.tnode,2]*2] <- 999
    
    prop.rules <- rules
    prop.rules[[t.node[prop.tnode,1]]][t.node[prop.tnode,2]] <- prop.rule
    
    ### Prune step
    temp.node <- c(0, 0)
    for(i in 2:length(prop.i.nodes)){
        temp.t <- grep(999, prop.i.nodes[[i]])
        if(length(temp.t) > 0){
            temp.node <- rbind(temp.node, cbind(i, temp.t))
        }
    }
    temp.node <- temp.node[-1,]
  
    # find internal nodes with two terminal child nodes
    count <- 0
    for(i in 2:dim(temp.node)[1]){
        if(temp.node[i,2] %% 2 == 0){if(temp.node[i-1,2] == (temp.node[i,2]-1) & temp.node[i,1] == temp.node[(i-1),1]){count <- count + 1}}
    }
    
    # Transition ratio
    TRANS <- log(p.prune) + log(num.t.node) - log(max(prop.prob[prop.pred]/sum(prop.prob[enough.unique]),0)) + log((unique.len-1)) - log(p.grow) - log(count)
  
    # Likelihood ratio
    nlL <- length(R.L)
    nlR <- length(R.R)
  
    LH <- log(sqrt(sigma2*(sigma2+(nlL+nlR)*sigma_mu))/sqrt((sigma2+nlL*sigma_mu)*(sigma2+nlR*sigma_mu))) + (sigma_mu / (2*sigma2) * ( (sum(R[R.L]))^2/(sigma2+nlL*sigma_mu) +(sum(R[R.R]))^2/(sigma2+nlR*sigma_mu) - (sum(R[union(R.R, R.L)]))^2/(sigma2+(nlR+nlL)*sigma_mu) ) )
  
    # Structure ratio
    d <- t.node[prop.tnode,1]-1
    STR <- log(alpha) + 2*log((1- alpha / (2 + d)^beta )) - log((1 + d)^beta - alpha) + log(max(prop.prob[prop.pred]/sum(prop.prob[enough.unique]),0)) - log((unique.len-1))
  
    r <- TRANS+LH+STR
    if(r > log(runif(1))){
        i.nodes <- prop.i.nodes;
        rules <- prop.rules;
    }
    
    t.node <- c(0, 0)
    for(i in 2:length(i.nodes)){
        temp <- grep(999, i.nodes[[i]])
        if(length(temp) > 0){
            t.node <- rbind(t.node, cbind(i, temp))
        }
    }
    t.node <- t.node[-1,]  # Level and Location of each terminal node
    
    # add the next level to i.nodes and rules
    if(t.node[dim(t.node)[1],1]==length(i.nodes)){
        i.nodes[[length(i.nodes)+1]] <- rep(NA, 2^(length(i.nodes)));
        rules[[length(rules)+1]] <- rep(NA, 2^(length(rules)))
    }
    
    
  return(list(i.nodes=i.nodes, rules=rules))
}
