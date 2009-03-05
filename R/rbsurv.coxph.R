################################################################################
#Computing negative loglik and AIC after fitting a Cox PH model
#When there are ties, use Efron (1977), Brelow (1974), or Cox (1972)
rbsurv.coxph.ties <- function(x.train, y.train, st.train, x.test, y.test, st.test, method="efron", k=2){
     

     if(!is.matrix(x.train)) x.train <- as.matrix(x.train)
     #Check errors
     tmp <- suppressWarnings(class(try(coxph(Surv(y.train, st.train) ~ x.train,  method=method), silent=TRUE)))
     if(tmp=="try-error")   return(list(nloglik = NA, AIC = NA, BIC = NA))

     fit <- suppressWarnings(coxph(Surv(y.train, st.train) ~ x.train,  method=method))
     b <- fit$coefficients
     #print(fit)

     i <- sort.list(y.test)
     y.test <- y.test[i]
     st.test <- st.test[i]
     x.test <- x.test[i,]

     if(!is.matrix(x.test))  x.test  <- as.matrix(x.test)
     p <- ncol(x.test)      # N of genes
     n <- length(y.test)    # N of samples
     tD <- unique(y.test[which(st.test==1)]) # distinct event times
     D <- length(tD)       # N of distinct event times

     #Compute log-likelihood
     si <- matrix(NA, D, p)
     di <- c()
     for(i in 1:D) {
        kk <- which((st.test==1) & (y.test==tD[i]))
        di <- c(di, length(kk))
        if(length(kk) >1) {
           if(is.null(dim((x.test[kk,])))) si[i,] <- sum(x.test[kk,])
           else si[i,] <- apply(x.test[kk,], 2, sum)   
        }  
        if(length(kk)==1) si[i,] <- x.test[kk,]    
     }

     if(method=="breslow") {
        tmp <- exp(x.test%*%b)
        LL <- c()
        for(i in 1:D) LL <- c(LL, sum(tmp[which(y.test >= tD[i])]))  
        LL <- sum(si%*%b) - sum(di*log(LL))
     }

     if(method=="exact") {
        tmp <- exp(x.test%*%b)
        sstar <- c()
        for(i in 1:D) sstar <- c(sstar, sum(tmp[which(y.test >= tD[i])]))  
        LL <- sum(si%*%b) - sum(log(sstar))
     }

     if(method=="efron") {
        tmp <- exp(x.test%*%b)
        tmp1 <- c(); tmp2 <- c()
        for(i in 1:D) {
            tmp1 <- c(tmp1, sum(tmp[which(y.test >= tD[i])])) 
            tmp2 <- c(tmp2, sum(tmp[which((st.test==1) & (y.test==tD[i]))])) 
        }
         
        LL <- sum(si%*%b)
        for(i in 1:D) for(j in 1:di[i]) LL <- LL - log(tmp1[i]-(j-1)/di[i]*tmp2[i]) 
     }


     #Compute AIC
     if(sum(is.na(b)) >0) {
        LL <- NA
        AIC <- NA
        BIC <- NA
      }
      else
      { 
        AIC <- -2*LL+k*p 
        BIC <- -2*LL+p*log(n) 
      }

     return(list(nloglik=-LL, AIC=AIC, BIC=BIC))
}

################################################################################
#Computing Loglik and AIC after fitting a Cox PH model
#When there no ties
rbsurv.coxph <- function(x.train, y.train, st.train, x.test, y.test, st.test, k=2){

     if(!is.matrix(x.train)) x.train <- as.matrix(x.train)

     #Check errors
     tmp <- suppressWarnings(class(try(coxph(Surv(y.train, st.train) ~ x.train), silent=TRUE)))
     if(tmp=="try-error")   return(list(nloglik = NA, AIC = NA, BIC=NA))

     fit <- suppressWarnings(coxph(Surv(y.train, st.train) ~ x.train))
     b <- fit$coefficients
     #print(fit)

     i <- sort.list(y.test)
     y.test <- y.test[i]; st.test <- st.test[i]; x.test <- x.test[i,]

     if(!is.matrix(x.test))  x.test  <- as.matrix(x.test)
     p <- ncol(x.test)      # N of genes
     n <- length(y.test)    # N of samples
     D <- which(st.test==1) # distinct event times ###NO D???
     n.D <- length(D)       # N of distinct event times

     #Compute log-likelihood
     s <- (x.test %*% b)
     exps <- exp(s)
     LL  <- sum(s[D])
     for(i in 1:n.D) LL <- LL-log(sum(exps[D[i]:n]))

     #Compute AIC
     if(sum(is.na(b)) >0) {
        LL <- NA
        AIC <- NA
        BIC <- NA
      }
      else
      { 
        AIC <- -2*LL+k*p 
        BIC <- -2*LL+p*log(n) 
      }

     return(list(nloglik=-LL, AIC=AIC, BIC=BIC))
}

################################################################################
