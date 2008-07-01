##########################################################################
#
#        
#    Robust Likelihood-based Survival Modeling for High-throughput Data
#
#             by
#
#      HyungJun Cho, Sukwoo Kim, and Soo-heang Eo
#      Deparments of Statistics and Biostatistics
#      Korea University
#
#        June 2008
#
##########################################################################

.First.lib <- function(lib, pkg) {  
     invisible()
     if(.Platform$OS.type=="windows" && require(Biobase) && interactive() 
        && .Platform$GUI=="Rgui") { addVigs2WinMenu("rbsurv") }
}

###################################################################################
#
# Main function  for users         
#
###################################################################################
rbsurv <- function(time, status, x, z=NULL, alpha=1, gene.ID=NULL, method="efron", 
                              n.iter=10, n.fold=3,  n.seq=1, seed = 1234, max.n.genes=nrow(x))
{



     ##########
     #Preparation
     set.seed(seed)
     n.samples <- ncol(x) 
     n.genes    <- nrow(x) 
     if(is.data.frame(x)==FALSE) x <- data.frame(x)
     colnames(x) <- 1:ncol(x)           
     rownames(x) <- 1:nrow(x) 
     x <- t(x) #convert (gene x sample) into (sample x gene)

     if(is.null(gene.ID)==TRUE)  gene.ID <-  as.character(1:n.genes)
     nfold <- max(2, min(n.fold, n.samples)) # 2 ~ N

     if((n.genes < 5) | (n.samples < 10)) stop('Too few genes or samples')
     require(survival)
     cat("Please wait...")

     #############
     #Reduce genes
     gene.ID.sub  <- gene.ID
     if(n.genes > max.n.genes)
    {
           tstat <- c() 
           for(i in 1:n.genes) {
                fit <- coxph(Surv(time, status) ~ x[,i],  method=method)
                tstat <- c(tstat, summary(fit)$logtest[1])
           }

           k <- which(length(tstat)-rank(tstat)+1 <= max.n.genes) 
           x <- x[,k] 
           gene.ID.sub  <- gene.ID[k]
     }

     #############
     #Select significant covariates
      covariates <- " NONE"
      if((is.null(z)==FALSE) | (alpha < 1)) 
      {
          k.covariates <- sig.covariates(time=time, status=status, z=z, method=method, alpha=alpha)
          if(length(k.covariates) >0) {
               z <- z[ , k.covariates, drop=FALSE]
               covariates <- colnames(z) 
          }
          if(length(k.covariates) ==0)  z <- NULL
       }

     ###############
     #Survival modeling
     gene.list <- c()
     probe.ID.sub <- probe.ID <- as.character(1:ncol(x))
     xx <- x
     out.model <- data.frame(matrix(NA, 1,  6))
     colnames(out.model) <- c("Seq","Order","Gene","nloglik","AIC",  "Selected")

     for(i in 1:n.seq) 
    {
            #print("n.seq"); print(i)
            if(ncol(x) < 5) break

            out <- rbsurv.sub(time=time, status=status, x=xx,  z=z,  
                                        gene.ID=probe.ID.sub, 
                                        method=method, 
                                        n.iter=n.iter, nfold=nfold)
            Seq <- rep(i, nrow(out$model))
            out.model <- rbind(out.model, cbind(Seq, out$model)) 

            k <- which.min(out$model$AIC) 
            if(k >1) gene.list <- c(gene.list, out$model$Gene[2:k]) 

            xx <- x
            if(length(gene.list) > 0) {
               k <- match(gene.list, probe.ID)
               xx <- x[,-k]
               probe.ID.sub <- probe.ID[-k]
            }
     }
     out.model <- out.model[-1,]

     ##########
     #ID matching
     ip <- which(as.numeric(out.model$Gene)>0)
     if(length(ip) >0) {
        k <- as.numeric(out.model$Gene[ip])
        out.model$Gene[ip] <- gene.ID.sub[k]
     }
     gene.list <- gene.ID.sub[as.numeric(gene.list)]

     cat(" Done. \n")
     return(list(n.genes=n.genes, n.samples=n.samples, 
                     method=method, 
                     n.iter=n.iter, n.fold=n.fold,  covariates=covariates,
                     model=out.model, gene.list=gene.list)) 

}


######################################################################
#Select signficant covariates
sig.covariates <- function(time=time, status=status, z=z, method=method, alpha=0.05)
{
        pvalue <- c() 
        for(i in 1:ncol(z)) {
             fit <- coxph(Surv(time, status) ~ z[,i],  method=method)
             pvalue <- c(pvalue, summary(fit)$logtest[3])
       }
       k <- which(pvalue < alpha) 
       return(k)
}


#END################################################################
