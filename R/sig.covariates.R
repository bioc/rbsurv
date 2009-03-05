'sig.covariates' <- function(time=time, status=status, z=z, method=method, alpha=0.05)
{
       pvalue <- c()
       for(i in 1:ncol(z)) {
             fit <- coxph(Surv(time, status) ~ z[,i],  method=method)
             pvalue <- c(pvalue, summary(fit)$logtest[3])
       }
       k <- which(pvalue < alpha)
       return(k)
}
