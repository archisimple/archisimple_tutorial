

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
source("src/io_archisimple.R")

### 1 . Read the parameter file

params <- read_archisimple("src/archisimple93/paramarch93.txt")

### 2 . Modify parameters

params$sim_length <- 10 ## put a minimal value of 10 here


### 3 . Save the parameters back to a text file to be read by ArchiSimple

write_archisimple(params, "src/archisimple93/paramarch93.txt")


### 4. run archisimple

setwd("src/archisimple93/")
system("./archisimple")
setwd("../../")


### 5. Get the output from the simulation

rs <- fread("src/archisimple93/seg.txt", header = T)

### 6. Plot the output

## From the side
rs %>%
  ggplot() +
  theme_classic() +
  geom_segment(aes(x = X1, y = -Z1, xend = X2, yend = -Z2), alpha=0.9) +
  coord_fixed()

# From the top
rs %>%
  ggplot() +
  theme_classic() +
  geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2), alpha=0.9) +
  coord_fixed()


