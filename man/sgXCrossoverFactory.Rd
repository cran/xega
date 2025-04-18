% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sgXCrossoverFactory.R
\name{sgXCrossoverFactory}
\alias{sgXCrossoverFactory}
\title{Factory for configuring a gene-dependent Crossover function.}
\usage{
sgXCrossoverFactory(algorithm = "sga", method = "CrossGene")
}
\arguments{
\item{algorithm}{Specifies algorithm. 
Available: "sga", "sgde", "sgperm", "sge", sgp". 
Default: "sga".}

\item{method}{Crossover method.  Algorithm (gene representation) 
dependent. Default: \code{CrossGene()}. 
Must be available in the gene-specific 
crossover factories.}
}
\value{
Crossover function from the crossover factory of 
                   the selected package.
}
\description{
\code{sgXCrossoverFactory()} selects 
             \enumerate{
             \item the algorithm-specific crossover factory and 
             \item the method in this factory. 
             }
}
\details{
The available methods for each algorithm are: 
   \itemize{
    \item "sga": 
       "Cross2Gene", "UCross2Gene", "UPCross2Gene",
       "CrossGene", "UCrossGene", "UPCrossGene".
    \item "sge": 
       "Cross2Gene", "UCross2Gene", "UPCross2Gene",
       "CrossGene", "UCrossGene", "UPCrossGene".
    \item "sgede": 
       "CrossGene", "UCrossGene", "UPCrossGene".
    \item "sgp": 
       "CrossGene", "Cross2Gene", "AllCrossGene", "AllCross2Gene",
       "FilterCrossGene", "FilterCross2Gene".
    \item "sgde": 
       "CrossGene", "UCrossGene", "UPCrossGene".
    \item "sgperm": 
       "CrossGene", "Cross2Gene".
 }
}
\examples{
sgXCrossoverFactory(algorithm="sga", method="CrossGene")

}
\seealso{
Other Configuration: 
\code{\link{sgXDecodeGeneFactory}()},
\code{\link{sgXGeneMapFactory}()},
\code{\link{sgXInitGeneFactory}()},
\code{\link{sgXMutationFactory}()},
\code{\link{sgXReplicationFactory}()}
}
\concept{Configuration}
