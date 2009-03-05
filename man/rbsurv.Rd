\name{rbsurv}
\alias{rbsurv}
\title{Robust likelihood-based survival modeling}
\description{
This selects survival-associated genes with microarray data.
}
\usage{
        rbsurv(time, ...)
}

\arguments{

\item{time}{an object for which the extraction of model rbsurv is meaningful.}
\item{...}{other arguments}
}


\references{
Cho,H., Yu,A., Kim,S., Kang,J., and Hong S-M. (2009). 
Robust likelihood-based survival modeling for microarray gene expression Data, 
\emph{Journal of Statistical Software}, 29(1):1-16. URL http://www.jstatsoft.org/v29/i01/.
 
 }

\author{
   HyungJun Cho, Sukwoo Kim, Soo-heang Eo, and Jaewoo Kang
}

\seealso{
    \code{\link{rbsurv.default}}
}

\examples{

library(rbsurv)
data(gliomaSet)
x <- exprs(gliomaSet)
x <- log2(x)
time <- gliomaSet$Time
status <- gliomaSet$Status
z <- cbind(gliomaSet$Age, gliomaSet$Gender) 

fit <- rbsurv(time=time, status=status, x=x,  method="efron", max.n.genes=20, n.iter=10, n.fold=3,  n.seq=1)
fit$model

}

\keyword{models}


