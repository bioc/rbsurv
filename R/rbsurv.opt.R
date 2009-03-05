#Choose a rule
option.rbsurv <- function(x.train, y.train, st.train, x.test, y.test, st.test, tie, method, k=2){

     if(tie==FALSE) tmp <- rbsurv.coxph(x.train, y.train, st.train, x.test,  y.test, st.test, k)
     if(tie==TRUE)  tmp <- rbsurv.coxph.ties(x.train, y.train, st.train, x.test,  y.test, st.test, method, k)
     return(list(nloglik=tmp$nloglik, AIC=tmp$AIC, BIC=tmp$BIC))
}


option.rbsurv.null <- function(y.train, st.train, y.test, st.test, tie, method){

     if(tie==FALSE) tmp <- rbsurv.coxph.null(y.train, st.train, y.test, st.test)
     if(tie==TRUE)  tmp <- rbsurv.coxph.ties.null(y.train, st.train, y.test, st.test, method)
     return(list(nloglik=tmp$nloglik, AIC=tmp$AIC, BIC=tmp$BIC))
}
