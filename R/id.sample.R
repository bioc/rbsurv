'id.sample' <- function(x, nfold=3){
     n1 <- sum(x)
     n0 <- sum(1-x)

     if(n1 < 5) stop(' The number of uncensored cases is too small!')
     if(n1 < 3*nfold) stop(' The number of fold is too large!')
     if(n0 <=0) {
        id <- sample(rep((1:nfold),n1)[1:n1])
        return(id)
     }

     id1 <- sample(rep((1:nfold),n1)[1:n1])
     id0 <- sample(rep((1:nfold),n0)[1:n0])
     id <- c(id0, id1)

     i <- (1:length(id))[sort.list(x)]
     id <- id[sort.list(i)]
     #print(table(id, x))
     return(id)

}
