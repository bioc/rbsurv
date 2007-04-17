################################################################################
rbsurv.coxph.ties.null <- function(y.train, st.train, y.test, st.test, method="efron"){
     
     i <- sort.list(y.test)
     y.test <- y.test[i]
     st.test <- st.test[i]

     n <- length(y.test)    # N of samples
     tD <- unique(y.test[which(st.test==1)]) # distinct event times
     D <- length(tD)       # N of distinct event times

     #Compute log-likelihood
     si <- 0
     di <- c()
     for(i in 1:D) {
        kk <- which((st.test==1) & (y.test==tD[i]))
        di <- c(di, length(kk))  
     }

     if(method=="breslow") {
        tmp <- rep(1, n)
        LL <- c()
        for(i in 1:D) LL <- c(LL, sum(tmp[which(y.test >= tD[i])]))  
        LL <- -sum(di*log(LL))
     }

     if(method=="exact") {
        tmp <- rep(1, n)
        sstar <- c()
        for(i in 1:D) sstar <- c(sstar, sum(tmp[which(y.test >= tD[i])]))  
        LL <- -sum(log(sstar))
     }

     if(method=="efron") {
        tmp <- rep(1, n)
        tmp1 <- c(); tmp2 <- c()
        for(i in 1:D) {
            tmp1 <- c(tmp1, sum(tmp[which(y.test >= tD[i])])) 
            tmp2 <- c(tmp2, sum(tmp[which((st.test==1) & (y.test==tD[i]))])) 
        }
         
        LL <- 0
        for(i in 1:D) for(j in 1:di[i]) LL <- LL - log(tmp1[i]-(j-1)/di[i]*tmp2[i]) 
     }

     #Compute AIC
     AIC <- -2*LL
     BIC <- -2*LL

     return(list(nloglik=-LL, AIC=AIC, BIC=BIC))
}

##########################################################################
rbsurv.coxph.null <- function(y.train, st.train, y.test, st.test){

     i <- sort.list(y.test)
     y.test <- y.test[i]
     st.test <- st.test[i]

     n <- length(y.test)    # N of samples
     D <- which(st.test==1) # distinct event times ###NO D???
     n.D <- length(D)       # N of distinct event times

     #Compute log-likelihood
     exps <- rep(1, n)
     LL  <- 0
     for(i in 1:n.D) LL <- LL-log(sum(exps[D[i]:n]))
     AIC <- -2*LL
     BIC <- -2*LL

     return(list(nloglik=-LL, AIC=AIC, BIC=BIC))
}

##########################################################################
