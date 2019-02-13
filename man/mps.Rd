\name{mps}
\alias{mps}
\title{Calculate minimum patch size}
\encoding{UTF-8}
\description{Calculates the minimum patch size below which a population cannot persist}
\usage{mps(r1 = 0.4631507, sigma2 = 285.6852, alpha = 1)}
\arguments{
  \item{r1}{rate of population growth in favourable habitat patches}
  \item{sigma2}{dispersal variance (of the dispersal distance of an individual)}
  \item{alpha}{ratio of the largest to the smallest side of a rectangular habitat patch (in the case of a square, alpha = 1)}
}
\details{Calculates the minimum size of a favourable habitat patch delimited by linear infrastructures such that the population size remains constant or increases}
\value{
  \item{Lmin}{The minimum size of the smallest side of a rectangular patch (i.e., the smallest distance between linear infrastructures) below which a population cannot persist}
  \item{Amin}{The minimum area of a patch below which a population cannot persist}
}
\references{
  Borda-de-Água L, Navarro L, Gavinhos C, Pereira HM. 2011. Landscape Ecology 26: 253-265. 
  
  Ceia-Hasse A, Borda-de-Água L, Grilo C, Pereira HM. 2017. Global Ecology and Biogeography 26: 592-600.
}
\author{Ana Ceia-Hasse, Luís Borda-de-Água, Henrique M. Pereira}
\examples{
# Run mps with default parameter values:

  mps()
  # mps(r1 = 0.4631507, sigma2 = 285.6852, alpha = 1)
  # Lmin =  78.02479 km
  # Amin =  6087.867 km2
  # r1 and sigma2 values for Puma concolor
}