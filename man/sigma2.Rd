\name{sigma2}
\alias{sigma2}
\title{Calculate dispersal variance}
\encoding{UTF-8}
\description{Calculates the dispersal variance of the dispersal distance of an individual}
\usage{sigma2(sigma_median = 7*sqrt(129.89), u = 1/16)}
\arguments{
  \item{sigma_median}{dispersal median (km/generation)}
  \item{u}{mortality rate (per year)}
}
\details{
  "u" was obtained from mean life span (mls) (u = 1/mls)
  
  "sigma_median" was obtained from home range (HR; km2) (sigma_median = 7*sqrt(HR) (mammals; Bowman et al. 2002. Ecology 83: 2049-2055); sigma_median = 12*sqrt(HR) (birds; Bowman 2003. Canadian Journal of Zoology 81: 195-202))
}
\value{dispersal variance (km2/year)}
\references{
  Pereira HM, Daily GC. 2006. Ecology 87: 1877-1885. 
  
  Ceia-Hasse A, Borda-de-Água L, Grilo C, Pereira HM. 2017. Global Ecology and Biogeography 26: 592-600.
}
\author{Ana Ceia-Hasse, Luís Borda-de-Água, Henrique M. Pereira}
\examples{
# Run sigma2 with default parameter values:

  sigma2()
  # sigma2(sigma_median = 7*sqrt(129.89), u = 1/16)
   # default values for Puma concolor (Jones et al. 2009. Ecology 90: 2648, Pereira and Daily 2006)
  # [1] 285.6852
}