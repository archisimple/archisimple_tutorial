

# The code should compile with any c++11 compiler, e.g. for g++: MinGW has been tested. Then create a new system environment variable. Path --> C:\MinGW\bin
# 
# Open the terminal:
#   cd src/archisimple93/
# 
# FOR MAC AND LINUX
#   g++ *.cpp -std=c++11 -o archisimple   
#   archisimple   
#
#
# FOR WINDOWS
#   g++ *.cpp -std=c++11 -o archisimple.exe   
#   archisimple.exe

library(tidyverse)
library(data.table)
setwd("~/Dropbox/science/projects/archisimple/archisimple_tutorial/2_batch_exemple//")
source("src/io_archisimple.R")

### 1 . Read the parameter file

params <- read_archisimple_xml("base_parameter.xml")


### 2 : create parameter ranges to run the simulations in batch

P_sim_length <- c(10,20)
P_dmax <- seq(1.2, 2, 0.1)


### 3 . Loop over the parameter space

all_sims <- NULL # Table to containn all the simulation data
inc <- 0 # increment for the simulations
setwd("src/archisimple93/")
for(dmax in P_dmax){
  for(sim_length in P_sim_length){
    
    ## increment 
    inc <- inc+1
    print(paste0("simulation ",inc))
    
    ### 4. Update the parameter in the loop
    params$sim_length <- sim_length ## put a minimal value of 10 here
    params$dmax <- dmax ## put a minimal value of 10 here
    
    ### 3 . Save the parameters back to a text file to be read by ArchiSimple
    write_archisimple_XML(params, path = "parameter.xml")
    
    ### 4. run archisimple
    system("./archisimple")
    
    ### 5. Get the output from the simulation
    rs <- fread("myroot.txt", header = T) %>% 
      mutate(sim_length = sim_length) %>%
      mutate(dmax = dmax) %>%
      mutate(sim = inc) # add a column with the simulation id
    
    ### 6. Store all the input in one big table
    all_sims <- rbind(all_sims, rs)
  }
}
setwd("../../")  



### 6. Plot the output of all the simulations (/!\ takes time for large simulations)

## From the side
all_sims %>%
  ggplot() +
  theme_classic() +
  geom_segment(aes(x = X1, y = -Z1, xend = X2, yend = -Z2), alpha=0.9) +
  coord_fixed() + 
  facet_grid(dmax~sim_length)

# From the top
rs %>%
  ggplot() +
  theme_classic() +
  geom_segment(aes(x = X1, y = -Z1, xend = X2, yend = -Z2), alpha=0.9) +
  coord_fixed()


