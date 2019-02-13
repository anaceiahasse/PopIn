\name{indmodel}
\alias{indmodel}
\title{Simulate source-sink population dynamics using a spatially-explicit individual-based model}
\encoding{UTF-8}
\description{Simulates population dynamics in landscapes with different types of habitat that have different habitat qualities (sources and sinks), including landscapes fragmented by linear infrastructures, using a spatially-explicit individual-based model that keeps track of the features of the landscape and of the population}
\usage{indmodel(land_r = "", param = "", file = "")}
\arguments{
\item{land_r}{A matrix that corresponds to the landscape}
\item{param}{A list with the model parameter values. See also "examples". The parameters are the following:
	
	1. n_rep: the number of replicates to be run
	
	2. hab_aff: the affinity of the species for each type of habitat in the landscape, in increasing order; e.g., if there are two types of habitat in the landscape matrix, "0" and "1" (corresponding to lower and higher quality habitats, respectively), and the affinity is 0.1 for the lower quality habitat ("0") and 0.9 for the higher quality habitat ("1"), then hab_aff = c(0.1, 0.9)
	
	3. n_steps: the number of steps in the simulation
	
	4. init_population: the population size at the start of the simulation
	
	5. hr_size: the size (number of cells) of the home range
	
	6. birth_rate: the maximum number of female juveniles produced by a breeding female (maximum fecundity)

	7. breeding_age: the age at which individuals start breeding

	8. survival: the probability of survival

	9. distance_weight: a higher distance weight corresponds to a higher contribution of cells closer to the home range centroid (thus the optimal home range shape is as close to a circle as possible)

	10. dispersal_distance: the distance (number of cells) that an individual has to disperse in order to settle a home range

	11. dispersal_mode: the dispersal mode; case 0: global dispersal; case 1: local dispersal, habitat search in a local kernel; case 2: local dispersal, random walk

	12. sink_avoidance: the probability that an individual avoids dispersing through cells that correspond to sink habitat

	13. neigh_avoidance: the probability that an individual avoids dispersing through cells that are already occupied

	14. sink_mortality: the probability that an individual dies while dispersing through cells that correspond to sink habitat

	15. plotYN: draws plots with results if set to TRUE
	
	16. writeYN: creates a directory containing files with results if set to TRUE

	17. dir_name: the name of the directory where files storing the results are saved

	18. file_name: the name of the file that will store, for each replicate: the values of the parameters, the landscape used, and for each time step: the position of each individual in the landscape (hrmaphist), the age of each individual (ageshist), and the population size (popsize)
}
\item{file}{A file with the model parameter values. See above and default file provided with the package (pim.txt)}
}
\details{
The function provides a default landscape and default parameter values using a file ("pim.txt"). It can read a list or a file with the parameter values, but if both a list and a file are given, an error message is returned.
}
\value{
\item{$popsize_av_end}{The population size (number of individuals) at the end of the simulation, averaged across all replicates}
\item{$prob_ext}{The probability of population extinction (the number of populations that went extinct divided by the total number of replicates)}
\item{$prob_persist}{The probability of population persistence (the number of populations that persisted divided by the total number of replicates (or 1-prob_ext))}
\item{$time_ext_av}{The average time to extinction (number of time steps it took for a population to go extinct, averaged across the replicates in which populations went extinct)}
\item{$popsizematrix}{A matrix with the population size at each time step of the simulation (columns), for each replicate (lines)}
\item{$popsize_av_time_step}{The population size (number of individuals) at each time step of the simulation, averaged across all replicates}
	
When plotYN = T, a plot is created showing the population size at each time step, averaged across all replicates.
}
\author{Ana Ceia-Hasse, Luís Borda-de-Água, Henrique M. Pereira}
\examples{
# Run indmodel with default parameter values:

  indmodel()    # indmodel(land_r = "", param = "", file = "")
  # uses default landscape and default file with parameter values
  
# Use with your own landscape and/or parameter values, e.g.:
  # indmodel(param = list(
  #    n_rep = 2, hab_aff = c(0.0, 1.0), n_steps = 200, init_population = 10,
  #    hr_size = 1, birth_rate = 2.0, breeding_age = 1, survival = 0.4, 
  #    distance_weight = 0.001, dispersal_distance = 5.0, dispersal_mode = 2,
  #    sink_avoidance = 0.1, neigh_avoidance = 1.0, sink_mortality = 0.9, 
  #    plotYN = "TRUE", writeYN = "TRUE", 
  #    dir_name = "Res_mylist", file_name = "indmodel_rep"))
}