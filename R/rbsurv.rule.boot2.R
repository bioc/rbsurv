##########################################################################
#
# Robust likelihood-based selection
#
##########################################################################
rbsurv.rule.boot2 <- function(z.train, x.train, y.train, st.train, z.test=NULL, x.test=NULL, y.test=NULL, st.test=NULL, 
                                            n.iter=10, nfold=3, method="efron") {
       

     ##########
     #Preparation
     n.gene <- ncol(x.train)
     n.sample.train <- nrow(x.train)
     colnames(x.train) <- 1:n.gene

      if(is.null(y.test)==FALSE)  
      { 
         n.sample.test  <- nrow(x.test)
         colnames(x.test)  <- 1:n.gene
      }

     ##################
     #Initiate
     opt.genes <- c()
     opt.nloglik    <- c()
     opt.AIC  <- c()
     opt.BIC  <- c()

     opt.nloglik.train    <- c()
     opt.AIC.train  <- c()
     opt.BIC.train  <- c()

     ##################
     #Pick 0-gene model
     xx.train <- NA
     xx.test  <- NA

     #Evaluate by train set
     tie <- TRUE; if(length(unique(y.train)) == length(y.train)) tie <- FALSE
     tmp <- option.rbsurv(z.train, y.train, st.train, z.train, y.train, st.train, tie=tie, method=method, k=2)
     null.nloglik.train <- tmp$nloglik
     null.AIC.train <- tmp$AIC
     null.BIC.train <- tmp$BIC

     #Evaluate by test set     
      if(is.null(y.test)==FALSE)  
      { 
          tie <- TRUE; if(length(unique(y.test)) == length(y.test)) tie <- FALSE
          tmp <- option.rbsurv(z.train, y.train, st.train, z.test, y.test, st.test, tie=tie, method=method, k=2)
          null.nloglik <- tmp$nloglik
          null.AIC  <- tmp$AIC
          null.BIC  <- tmp$BIC
     }


     ####################
     #Pick k-gene model (k >=1)
     i.stop <- 0 
     x.train.opt    <- z.train
     x.train.cand <- x.train 

     for(jj in 1:(n.gene-1)) {
        n.gene.cand <- n.gene-jj+1
        out <- matrix(0, n.iter, n.gene.cand)
        for(i in 1:n.iter) {
            id <- id.sample(st.train, nfold=nfold)
            y.tr <- y.train[id!=1]
            y.te <- y.train[id==1]
            st.tr <- st.train[id!=1]
            st.te <- st.train[id==1]
            tie <- TRUE; if(length(unique(y.te)) == length(y.te)) tie <- FALSE
            tmp <- -1########
            for(j in 1:n.gene.cand) {
                x.tr  <- data.frame(x.train.opt[id!=1, , drop=FALSE], x.train.cand[id!=1, j, drop=FALSE])
                x.te <- data.frame(x.train.opt[id==1, , drop=FALSE], x.train.cand[id==1, j, drop=FALSE])
                out[i,j] <- option.rbsurv(x.tr, y.tr, st.tr, x.te, y.te, st.tr, 
                                    tie=tie, method=method, k=2)$nloglik
                #print(j) 
                #print(out[i,j])
            }
        }

        ii <- which(apply(apply(out, 2, is.na), 2, sum) >0)
        if(length(ii) >= ncol(out)) break 
        if(length(ii) > 0) out[,ii] <- max(out, na.rm=TRUE)+1
        out.sum <- apply(out, 2, sum, na.rm = TRUE)
        if(length(unique(out.sum)) == 1) break

        pick.gene <- which.min(out.sum) 
        if(length(pick.gene) != 1) next #FIXED
        pick.gene <- as.numeric(colnames(x.train.cand)[pick.gene])
        opt.genes <- c(opt.genes, pick.gene)
        x.train.opt  <- cbind(z.train, x.train[, opt.genes,drop=FALSE])
        x.train.cand <- x.train[,-opt.genes,drop=FALSE]

        #Evaluate by train set
        xx.train <- cbind(z.train,  x.train[,opt.genes,drop=FALSE])
        tie <- TRUE; if(length(unique(y.train)) == length(y.train)) tie <- FALSE
        tmp <- option.rbsurv(xx.train, y.train, st.train, xx.train, y.train, st.train, 
                                         tie=tie, method=method, k=2)
        opt.nloglik.train <- c(opt.nloglik.train, tmp$nloglik)
        opt.AIC.train  <- c(opt.AIC.train, tmp$AIC)
        opt.BIC.train  <- c(opt.BIC.train, tmp$BIC)

        #Evaluate by test set
        if(is.null(y.test)==FALSE)  
       { 
            xx.test  <- cbind(z.test, x.test[,opt.genes,drop=FALSE])
            tie <- TRUE; if(length(unique(y.test)) == length(y.test)) tie <- FALSE
            tmp <- option.rbsurv(xx.train, y.train, st.train, xx.test,  y.test, st.test, 
                                             tie=tie, method=method, k=2)
            opt.nloglik <- c(opt.nloglik, tmp$nloglik)
            opt.AIC  <- c(opt.AIC, tmp$AIC)
            opt.BIC  <- c(opt.BIC, tmp$BIC)
        }

        #stop if min number of classes is less than 2
        if(sum(st.train)*(1-1/nfold) <= length(opt.genes)+1) break    

    }

    ####### 
    #Outputs
    if(is.null(y.test)==TRUE)  
    { 
        opt.nloglik  <- opt.nloglik.train 
        opt.AIC       <- opt.AIC.train 
        null.nloglik <-  null.nloglik.train 
        null.AIC      <- null.AIC.train 
    }

    i <- 1:length(opt.genes)
    if(length(opt.genes) >0) final.out <- data.frame(i, opt.genes, opt.nloglik.train, opt.AIC.train, opt.nloglik, opt.AIC)
    if(length(opt.genes) ==0) final.out <- data.frame(rep(0,6))
    colnames(final.out) <- c("Order","Gene","Train.nloglik","Train.AIC", "nloglik","AIC")

    i <- 0
    final.out.null <- data.frame(i, i, null.nloglik.train, null.AIC.train,  null.nloglik, null.AIC)
    colnames(final.out.null) <- c("Order","Gene","Train.nloglik","Train.AIC", "nloglik","AIC")
    final.out <- rbind(final.out.null, final.out)

    return(final.out)

}

#END######################################################