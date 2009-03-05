##########################################################################
#
# Sub-function
#
##########################################################################
rbsurv.sub <- function(time, status, x, z=NULL, gene.ID=NULL, method="efron", n.iter=10, nfold=3)

{

        ##########
        #Preparation
        if((ncol(x) < 5) | (nrow(x) < 10))
        {
             print("Too few genes or samples")
             return(list(model=NULL, gene.list=NULL))  
        }
         
        ########
        #Modeling
        if(is.null(z)==TRUE) {
           out <- rbsurv.rule.boot(x.train=x, y.train=time, st.train=status, 
                                               n.iter=n.iter, nfold=nfold, method=method)
        }
        if(is.null(z)==FALSE) { #with covariates
           out <- rbsurv.rule.boot2(z.train=z, x.train=x, y.train=time, st.train=status,
                                               n.iter=n.iter, nfold=nfold, method=method)
        }

        ##########
        #ID matching
        n <- length(out[,2]) 
        if(n >1) out[2:n, 2] <- gene.ID[out[2:n,2]]


        ########
        #Selection
        Select <- rep(" ", nrow(out))
        if(which.min(out$AIC) >1) 
        {
             Select[2:which.min(out$AIC)] <- "*       "
             gene.list <- out$Gene[2:which.min(out$AIC)] 
        }
        else
        { 
             Select[1] <- "*      "
             gene.list <- NULL
        }

        
        #######
        #Outputs
        out <- cbind(out[,1:4], Select)
        rownames(out) <- c(0:(nrow(out)-1))
        colnames(out) <- c("Order","Gene","nloglik","AIC",  "Selected")
        out$nloglik <- round(out$nloglik, 2)
        out$AIC  <- round(out$AIC, 2)
        return(list(model=out, gene.list=gene.list))  
}

##END#####################################################################

