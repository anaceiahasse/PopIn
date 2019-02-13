\name{mlid}
\alias{mlid}
\title{Calculate maximum linear infrastructure density}
\encoding{UTF-8}
\description{Calculates the maximum linear infrastructure density above which a population cannot persist}
\usage{mlid(r1 = 0.4631507, r0 = -(1/16)*10^3, li.width = 0.01)}
\arguments{
  \item{r1}{rate of population growth in favourable habitats}
  \item{r0}{rate of population decay in unfavourable habitats}
  \item{li.width}{the width of the linear infrastructure (km)}
}
\details{r0 corresponds to the rate of population decay due to direct mortality occurring when individuals cross the linear infrastructure (e.g., roadkills due to wildlife-vehicle collisions)}
\value{
  \item{Dmax}{The maximum linear infrastructure density above which a population cannot persist}
}
\references{
  Borda-de-Água L, Navarro L, Gavinhos C, Pereira HM. 2011. Landscape Ecology 26: 253-265.

  Ceia-Hasse A, Borda-de-Água L, Grilo C, Pereira HM. 2017. Global Ecology and Biogeography 26: 592-600.
}
\author{Ana Ceia-Hasse, Luís Borda-de-Água, Henrique M. Pereira}
\examples{
# Run mlid with default parameter values:

  mlid()
  # mlid(r1 = 0.4631507, r0 = -(1/16)*10^3, li.width = 0.01)
  # Dmax =  0.7355901 km/km2
  # r1 and r0 values for Puma concolor
}