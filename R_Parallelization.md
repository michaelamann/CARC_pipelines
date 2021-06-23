# Running R in Parallel 



In R, there are numerous packages that can be used to parallelize code (parallel, snow, foreach, etc.), but these packages use unique syntaxes and none of them work for all cases of parallelization. Instead, the [future package](https://github.com/HenrikBengtsson/future) which has a unified syntax can be used to utilize multiple processors and multiple nodes for iterative processes. Briefly (see [future package](https://github.com/HenrikBengtsson/future) for details), the future packagae split






## set up

``` bash

# create interactive node with 8 processors
qsub -I -lwalltime=4:00:00 -lnodes=1:ppn=8

# load R module
module load r-4.0.4-gcc-10.2.0-python3-dghog6f

```




## Future with for loop
``` R 
# Load the package 
library(future)

# custom function that waits a second, and then prints current step value
# will loop through this function as an example

slow_function <- function(step){
  
  # wait three seconds
  Sys.sleep(1)
  
  
  paste0("Step ", step, " completed")
}



#### running loop in serial ####
# timestampe before loop 
t1 <- proc.time()
# pre-allocate output 
output <- rep(NA, 50)

# for loop where it gets increasingly slower
for (i in 1:50){
  output[i] <- slow_function(i)

}
output

# timestampe after loop 
t2 <- proc.time()

# length of time
t2 - t1

 #  user   system  elapsed 
 #  0.038  0.005   50.098 



#### rerunning loop in parallel using all available cores ####
# rewrite loop so it works with future:
output <- list()

# going to send iterations to each of the 3 cores
plan(multicore)

t1 <- proc.time()
# for loop where it gets increasingly slower
for (i in 1:50){
  # writing each iteration code to y 
  
  output[[i]] <- future(
    # but code for each iteration with curly brackets
    {slow_function(i)}
  ) 
}

# evaluate, distribute each iteration here 
y <- value(y)
unlist(output)
# ouput is here


# timestampe after loop 
t2 <- proc.time()

# length of time
t2 - t1


 #  user  system  elapsed 
 #  1.694  1.661  12.014 
```





## Using future with tidyverse via the furrr package

``` R
library(future) # needed for parallelization
library(future.batchtools) # needed to submit jobs in parallel
library(furrr) # allows future to be used with tidyverse code
library(repurrrsive) # loads in data
library(tidyverse) 

# create a nested dataframe with each country as a row 
# we will iterate using the map function from the purrr package
country_nested <- 
  gap_simple %>%
  group_by(country) %>%
  nest() %>%
  ungroup() # this is necessary or furrr will be slow


custom_model <- 
  function(data){
    # just to slow it down so its more obvious it's in parallel
    Sys.sleep(.2)
    lm(lifeExp ~ pop + gdpPercap + year, data = data)
  }



# serial form (normal purrr)
t1 <- proc.time()
country_nested %>%
  mutate(lm_obj = map(data, custom_model))
t2 <- proc.time()
t2 - t1
 #  user  system  elapsed 
 #  0.547  0.042  29.058 


# using furrr package

# double check you have multiple cores available
availableCores()

# use all cores available
plan(multicore) 


t1 <- proc.time()
country_nested %>%
    mutate(lm_obj = future_map(data, custom_model)) # switched map to future_map!

# switch back to one core
plan(sequential)   
t2 <- proc.time()
t2 - t1

 #  user   system  elapsed 
 #  0.805  0.437   4.214 





````





## submitting jobs across muliple nodes 


``` R
library(Hmsc) # package for joint species distribution models
library(ape) # needed for phylogeny
library(future)
library(future.batchtools) # needed to submit pbs scripts


# bring in data for model
load("hmsc_setup.RData")


# Setting up the model
studyDesign = data.frame(Route = XData$Route)
rL = HmscRandomLevel(sData=xy)
XFormula = ~ hab + poly(clim,degree = 2,raw = TRUE)
TrFormula = ~Migration + LogMass

# Full model
m = Hmsc(Y=Y, XData = XData, XFormula=XFormula, 
         phyloTree = phyloTree, TrData = TrData, 
         TrFormula = TrFormula,
         distr="probit", studyDesign=studyDesign, 
         ranLevels=list(Route=rL))

# parameters
nChains = 4
nParallel = 4 
samples = 100



# set up that it will submit pbs scripts for each model
plan(sequential)

y <- list()
t1 <- proc.time()
# running 4 models, each with a different thinning value
for (thin in c(1,5,10,25)){
  y[[thin]] <- future({
  transient = 50*thin
  m = sampleMcmc(m, thin = thin, samples = samples, transient = transient,
                 nChains = 4, initPar = "fixed effects",
                 nParallel = nParallel)
                 
  # write model outputs to file               
  filename=file.path(paste0("Big_model_torque_chains_",as.character(nChains),"_samples_",as.character(samples),"_thin_",as.character(thin)))
  save(m,file=filename)
  }, seed = TRUE)
}
# evaluate expression
y <- value(y) 

t2 <- proc.time()
t2 - t1
# done!








# set up that it will submit pbs scripts for each model
plan(batchtools_torque)

y <- list()

# running 4 models, each with a different thinning value
for (thin in c(1,5,10,25)){
  y[[thin]] <- future({
  transient = 50*thin
  m = sampleMcmc(m, thin = thin, samples = samples, transient = transient,
                 nChains = 4, initPar = "fixed effects",
                 nParallel = nParallel)
                 
  # write model outputs to file               
  filename=file.path(paste0("Big_model_torque_chains_",as.character(nChains),"_samples_",as.character(samples),"_thin_",as.character(thin)))
  save(m,file=filename)
  }, seed = TRUE)
}
# evaluate expression
y <- value(y) 

```





