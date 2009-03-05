\name{rbsurv.default}
\alias{rbsurv.default}
\title{Robust likelihood-based survival modeling}
\description{
This selects survival-associated genes with microarray data.
}
\usage{
    \method{rbsurv}{default}(time, status, x, z=NULL, alpha=1, gene.ID=NULL, method="efron",
                   n.iter=10, n.fold=3,  n.seq=1,  seed=1234, max.n.genes=nrow(x),...)
}

\arguments{

\item{time}{a vector for survival times}
\item{status}{a vector for survival status, 0=censored, 1=event}
\item{x}{a matrix for expression values (genes in rows, samples in columns)}
\item{z}{a matrix for risk factors}
\item{alpha}{significance level for evaluating  risk factors; significant risk factors included with the alpha level 
    if alpha < 1}
\item{gene.ID}{a vector for gene IDs; if NULL, row numbers are assigned.}
\item{method}{a character string specifying the method for tie handling.  Choose one of "efron", "breslow",
        "exact". The default is "efron".  If there are no tied death times all the methods are equivalent. }
\item{n.iter}{the number of iterations for gene selection}
\item{n.fold}{the number of partitions of samples}
\item{n.seq}{the number of sequential runs or multiple models}
\item{seed}{a seed for sample partitioning}
\item{max.n.genes}{the maximum number of genes considered.
If the number of the input genes is greater than the given number, it is reduced by fitting individual Cox models.}
\item{...}{other arguments}

}

\value{
   \item{model}{ survival-associated gene model}
   \item{n.genes}{number of genes}
   \item{n.samples}{number of samples}
   \item{method}{method for tie handling}
   \item{covariates}{covariates}
   \item{n.iter}{number of iterations for gene seletion}
   \item{n.fold}{number of partitions of samples}
   \item{n.seq}{number of sequential runs or multiple models}
   \item{gene.list}{a list of genes included in the models}
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
\code{\link{rbsurv}}
}

\keyword{models}


