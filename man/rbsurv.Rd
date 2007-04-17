\name{rbsurv}
\alias{rbsurv}
\title{Robust likelihood-based survival modeling}
\description{
This selects survival-associated genes with microarray data.
}
\usage{
        rbsurv(time, status, x, z=NULL, alpha=0.05, gene.ID=NULL, method="efron", 
                   max.n.genes=100, n.iter=10, n.fold=3,  n.seq=1)
}

\arguments{

\item{time}{a vector for survival times}
\item{status}{a vector for survival status, 0=censored, 1=event}
\item{x}{a matrix for expression values (genes in rows, samples in columns)}
\item{z}{a matrix or vector for covariates }
\item{alpha}{significance level for covariates}
\item{gene.ID}{a vector for gene IDs; if NULL, row numbers are assigned.}
\item{method}{a character string specifying the method for tie handling.  Choose one of "efron", "breslow", 
        "exact". The default is "efron".  If there are no tied death times all the methods are equivalent. }
\item{max.n.genes}{The maximum number of genes considered. 
If the number of the input genes is greater than the given number, it is reduced by fitting individual Cox models. 
The default is 100.}
\item{n.iter}{number of iterations} 
\item{n.fold}{number of partitions of samples} 
\item{n.seq}{number of sequential runs}

}

\value{
   \item{model}{survival-associated gene model}
}

\references{
Cho H et al. Robust likelihood-based survival modeling for microarray gene expression data, submitted.
 }

\author{
   HyungJun Cho and Sukwoo Kim
}


\examples{

library(rbsurv)
data(glioma)
x <- log2(glioma.x) # do normalization if necessary
time <- glioma.y$Time
status <- glioma.y$Status

fit <- rbsurv(time=time, status=status, x=x,  method="efron", max.n.genes=20, n.iter=10, n.fold=3,  n.seq=1)
fit$model

}

\keyword{models}


