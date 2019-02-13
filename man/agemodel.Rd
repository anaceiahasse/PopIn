\name{agemodel}
\alias{agemodel}
\title{A spatially-implicit, age-structured, stochastic model to evaluate the impact of direct mortality caused by infrastructures on populations and their risk of extinction}
\encoding{UTF-8}
\description{Simulates the time evolution of a population with and without infrastructure driven mortality, enabling its study under different infrastructure mortality rates and carrying capacities, including its probability of extinction, and taking into account the specificities of the life history of the species. See "Details"}
\usage{agemodel(n.rep, param = "", file = "")}
\arguments{
\item{n.rep}{The number of replicates to be run}

\item{param}{A list with the model parameter values. See also "examples". The parameters are the following:
	
	1. years: the total number of years for which the simulation runs
	
	2. transient1: number of years; ensures that the age classes stabilize
	
	3. transient2: number of years; this period is used to estimate the average population size without infrastructure mortality (which is applied after transient2)
	
	4. transient3: number of years; ensures that the population "stabilizes" after infrastructure driven mortality starts being applied; only if the population survives until the end of transient3 then information is collected
	
	5. a: parameter of the Beverton-Holt relationship; used when n.rep = 1
	
	6. b: parameter of the Beverton-Holt relationship; used when n.rep = 1
	
	7. a.val: the values that "a" will take; this parameter is used when n.rep > 1

	8. f.mort: the fraction of the population killed by infrastructures per year; this parameter is used when n.rep = 1
	
	9. f.values: the values that the fraction of the population killed by infrastructures per year takes; this parameter is used when n.rep > 1
	
	10. n.ages: the number of age classes, i.e., the maximum number of years an individual can live
	
	11. min.age: the age at maturity, that is, the minimum age at which an animal can reproduce (in months)
	
	12. repr.period: the reproduction interval for a female (in months)
	
	13. months.of.recruitment: the months where animals reproduce
	
	14. prob.recr.month.m: in those months where the population can reproduce, the mean proportion that reproduces (logit scale)
	
	15. prob.recr.month.sd: in those months where the population can reproduce, the standard deviation of the proportion that reproduces (logit scale)
	
	16. prop.recr.age.m: the mean proportion per age that reproduces (linear scale)
	
	17. prop.recr.age.sd: the standard deviation of the proportion per age that reproduces (linear scale)
	
	18. clutch.size.m: the mean number of eggs per clutch (mean litter size; mean number of offspring)
	
	19. clutch.size.sd: the standard deviation of the number of eggs per clutch
		
	20. jsmm: the mean of juvenile survival rate for males
	
	21. jssdm: the standard deviation of juvenile survival rate for males
	
	22. jsmf: the mean of juvenile survival rate for females
	
	23. jssdf: the standard deviation of juvenile survival rate for females
	
	24. ysmm: the mean of yearling survival rate for males
	
	25. yssdm: the standard deviation of yearling survival rate for males
	
	26. ysmf: the mean of yearling survival rate for females
	
	27. yssdf: the standard deviation of yearling survival rate for females
	
	28. asmm: the mean of adult survival rate for males
	
	29. assdm: the standard deviation of adult survival rate for males
	
	30. asmf: the mean of adult survival rate for females
	
	31. assdf: the standard deviation of adult survival rate for females
	
	32. plotTF: draws plots with results if set to TRUE
	
	33. writeTF: creates a directory containing files with results if set to TRUE

	34. dirname.basic: the name of the directory where the files with results are stored when writeTF = TRUE and when n.rep = 1 
	
	35. dirname.rep: the name of the directory where the files with results are stored when writeTF = TRUE and n.rep > 1 
}
\item{file}{A file with the model parameter values. See above and default file provided with the package (pam.txt)}
}
\details{
The function takes as inputs the number of replicates (n.rep) and a list or a file with the remaining parameter values. If both a list and a file are given, an error message is returned. The package provides default parameter values using a file ("pam.txt").

Density dependence is introduced between the number of recruits and the number of eggs, thus the values of parameters a and b (Beverton-Holt relationship) should be chosen so that the number of recruits is always smaller than the number of eggs (i.e., a is smaller than or equal to b).
}
\value{
The function return the names of and the values used for the parameters. See "arguments" for the description of each parameter and "examples". When plotTF = T, plots are created; when writeTF = T, a directory is created which includes files with results (see below).

WHEN plotTF = T, FOUR PLOTS ARE CREATED:

	1. When n.rep = 1

		1.1. Plots 1 and 2: Time evolution of the simulated number of individuals during the years that the simulation runs, in January, and in December.

		1.2. Plot 3: Number of individuals killed by infrastructures per year from the time when such mortality is applied until the end of the simulation.

		1.3. Plots 4 and 5: Histogram of the logarithm of the simulated number of individuals killed per year during the infrastructure driven mortality regime, and a plot with the same data in double logarithmic scales.

	2. When n.rep > 1
	
		# the plots show the results of the simulations for different values of a.val, and several values of f.values:

		2.1. Plot a: Number of animals at the first month of recruitment.
			
		2.2. Plot b: Number of animals killed per year.

		2.3. Plot c: Probability of extinction.

		2.4. Plot d: Time to extinction (years).
			

WHEN writeTF = T, A DIRECTORY IS CREATED WHICH INCLUDES THE FILES:

	1. When n.rep = 1

		1.1. pop_size_jan: the time evolution of the simulated number of individuals in the population for the month of January, during the years that the simulation runs. The first column corresponds to the years and the second column corresponds to the number of individuals in the population in January, in each of the years.

  		1.2. pop_size_dec: the time evolution of the simulated number of individuals in the population for the month of December, during the years that the simulation runs. The first column corresponds to the years and the second column corresponds to the number of individuals in the population in December, in each of the years.

  		1.3. n_killed: the number of individuals killed by infrastructures per year from the time when such mortality is applied until the end of the simulation. The first column corresponds to the years and the second column corresponds to the number of individuals killed by infrastructures in each year.

	2. When n.rep > 1

  		2.1. a_values: values taken by a
  
  		2.2. f_values: values taken by the fraction of the population to be killed by infrastructures per year

  		2.3. population size before infrastructure mortality is applied (after the age classes have stabilized), for each value of a.val, averaged across replicates and f.values:
	
    		2.3.1. before recruitment (January):
	
      			2.3.1.1. n_annual_median: median population size
	
		  		2.3.1.2. n_annual_sd: standard deviation of the population size
	
		  		2.3.1.3. n_annual_mean: mean population size
	
		  		2.3.1.4. n_annual_max: maximum population size
	
		  		2.3.1.5. n_annual_min: minimum population size
	
    		2.3.2. after recruitment (December):
	
		  		2.3.2.1. n_annual_recr_median: median population size
	
		 		2.3.2.2. n_annual_recr_sd: standard deviation of the population size
	
		  		2.3.2.3. n_annual_recr_mean: mean population size
	
		  		2.3.2.4. n_annual_recr_max: maximum population size
	
		  		2.3.2.5. n_annual_recr_min: minimum population size

  		2.4. absolute maximum and minimum population size before infrastructure mortality is applied (after the age classes have stabilized), for each value of a.val:
	
    		2.4.1. before recruitment (January):
	
		  		2.4.1.1. n_annual_max_max: absolute maximum population size
	
		  		2.4.1.2. n_annual_min_min: absolute minimum population size
	
   			2.4.2. after recruitment (December):
	
		  		2.4.2.1. n_annual_recr_max_max: absolute maximum population size
	
		  		2.4.2.2. n_annual_recr_min_min: absolute minimum population size

  		2.5. population size after infrastructure mortality is applied, for the populations that persisted to the end of the simulation, for each combination of values of (a.val, f.values), averaged across replicates:
	
    		2.5.1. before recruitment (January):
	
		  		2.5.1.1. n_annual_m_median: median population size
	
		  		2.5.1.2. n_annual_m_sd: standard deviation of the population size
	
		  		2.5.1.3. n_annual_m_mean: mean population size
	
		  		2.5.1.4. n_annual_m_max: maximum population size
	
		  		2.5.1.5. n_annual_m_min: minimum population size
	
    		2.5.2. after recruitment (December):
	
		  		2.5.2.1. n_annual_recr_m_median: median population size
	
		  		2.5.2.2. n_annual_recr_m_sd: standard deviation of the population size
	
		  		2.5.2.3. n_annual_recr_m_mean: mean population size
	
		  		2.5.2.4. n_annual_recr_m_max: maximum population size
	
		  		2.5.2.5. n_annual_recr_m_min: minimum population size

  		2.6. absolute maximum and minimum population size after infrastructure mortality is applied, for the populations that persisted to the end of the simulation, for each combination of values of (a.val, f.values):
	
    		2.6.1. before recruitment (January):
	
		  		2.6.1.1. n_annual_m_max_max: absolute maximum population size
	
		  		2.6.1.2. n_annual_m_min_min: absolute minimum population size
	
    		2.6.2. after recruitment (December):
	
		  		2.6.2.1. n_annual_recr_m_max_max: maximum population size
	
		  		2.6.2.2. n_annual_recr_m_min_min: minimum population size

  		2.7. number of individuals killed by infrastructures for the populations that persisted to the end of the simulation, for each combination of values of (a.val, f.values), averaged across replicates:

    		2.7.1. n_killed_median: median number of individuals killed by infrastructures

    		2.7.2. n_killed_sd: standard deviation

    		2.7.3. n_killed_mean: mean

    		2.7.4. n_killed_max: maximum

    		2.7.5. n_killed_min: minimum

  		2.8. absolute maximum and minimum number of individuals killed by infrastructures for the populations that persisted to the end of the simulation, for each combination of values of (a.val, f.values):

    		2.8.1. n_killed_max_max: absolute maximum

    		2.8.2. n_killed_min_min: absolute minimum

  		2.9. f_mort_12: fraction of the population killed by infrastructures in December for the different combinations of values of (a.val, f.values)

  		2.10. results for the populations that went extinct, for the different combinations of values of (a.val, f.values)

    		2.10.1. prob_ext: probability of population extinction

    		2.10.2. mean_t: mean time to extinction

    		2.10.3. median_t: median time to extinction

    		2.10.4. sd_t: standard deviation of the time to extinction

    		2.10.5. max_t: maximum time to extinction

    		2.10.6. min_t: minimum time to extinction
}
\references{Borda-de-Água L, Grilo C, Pereira HM. 2014. Ecological Modelling 276: 29-37}
\author{Ana Ceia-Hasse, Luís Borda-de-Água, Henrique M. Pereira}
\examples{
# Run agemodel with default parameter values:
  # uses default file provided with the package
  
  agemodel(1,,)    # agemodel(n.rep = 1, param = "", file = "")
  
  agemodel(2,,)    # agemodel(n.rep = 2, param = "", file = "")
  
# Use with your own parameter values, e.g.:
  # agemodel(n.rep = 1, param = mylist_1rep)
  # agemodel(n.rep = 2, param = mylist_2rep)

  # mylist_1rep = list(
  #    years = 320, transient1 = 120, transient2 = 100, transient3 = 50, 
  #    f.mort = 0.05, a = 5, b = 5, 
  #    n.ages = 17, min.age = 37, repr.period = 18, 
  #    months.of.recruitment = c("May", "June", "July", "August"), 
  #    prob.recr.month.m = c(100,100,100,100,100,100,100,100,100,100,100,100), 
  #    prob.recr.month.sd = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
  #    prop.recr.age.m = rep(1, 17), prop.recr.age.sd = rep(0, 17), 
  #    clutch.size.m = 2.3, clutch.size.sd = 1, 
  #    jsmm = 0.86, jssdm = 0.09, jsmf = 0.86, jssdf = 0.06, 
  #    ysmm = 0.62, yssdm = 0.06, ysmf = 0.62, yssdf = 0.06, 
  #    asmm = 0.47, assdm = 0.05, asmf = 0.47, assdf = 0.05, 
  #    plotTF = T, writeTF = T, dirname.basic = "agemodel_1rep")

  # mylist_2rep = list(
  #    years = 320, transient1 = 120, transient2 = 100, transient3 = 50, 
  #    f.values = seq(0.1, 0.9, 0.1), a.val = seq(1, 5, 1), 
  #    n.ages = 17, min.age = 37, repr.period = 18, 
  #    months.of.recruitment = c("May", "June", "July", "August"), 
  #    prob.recr.month.m = c(100,100,100,100,100,100,100,100,100,100,100,100), 
  #    prob.recr.month.sd = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
  #    prop.recr.age.m = rep(1, 17), prop.recr.age.sd = rep(0, 17), 
  #    clutch.size.m = 2.3, clutch.size.sd = 1, 
  #    jsmm = 0.86, jssdm = 0.09, jsmf = 0.86, jssdf = 0.06, 
  #    ysmm = 0.62, yssdm = 0.06, ysmf = 0.62, yssdf = 0.06, 
  #    asmm = 0.47, assdm = 0.05, asmf = 0.47, assdf = 0.05, 
  #    plotTF = T, writeTF = T, dirname.rep = "agemodel_2rep")
}