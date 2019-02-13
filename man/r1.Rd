\name{r1}
\alias{r1}
\title{Calculate population growth rate in favourable habitats}
\encoding{UTF-8}
\description{Calculates the growth rate of a population in favourable habitats}
\usage{r1(b = 1.75, B = 2.75, s = (1 - (1/16)), maxage = 80)}
\arguments{
  \item{b}{0.5*litter size}
  \item{B}{breeding age}
  \item{s}{survival rate}
  \item{maxage}{maximum age}
}
\details{
  "survival rate" (s) can be computed from mortality rate (u) (s = 1 - u) or from mean life span (mls) (s = 1 - (1/mls))
  
  "maxage" can be computed from mean life span (mls) (maxage = 5 * mls)
}
\value{population growth rate in favourable habitats}
\references{
  Pereira HM, Daily GC, Roughgarden J. 2004. Ecological Applications 14: 730-742.
  
  Pereira HM, Daily GC. 2006. Ecology 87: 1877-1885. 
  
  Ceia-Hasse A, Borda-de-Água L, Grilo C, Pereira HM. 2017. Global Ecology and Biogeography 26: 592-600.
}
\author{Ana Ceia-Hasse, Luís Borda-de-Água, Henrique M. Pereira}
\examples{
# Run r1 with default parameter values:

  r1()
  # r1(b = 1.75, B = 2.75, s = (1 - (1/16)), maxage = 80) 
  # default values for Puma concolor (Pereira and Daily 2006)
  # [1] 0.4631507
}