
## IO functions to read and write archisimple parameter files
library(tidyverse)

### READ function

path <- "~/Dropbox/science/projects/archisimple/archishiny/www/paramarch93.txt"
options(scipen = 999) 
read_archisimple <- function(path){
  
  # names of the parameters
  var_names <- c("sim_length", 
                 "emission_rate_sem", 
                 "dmax_prop_sem", 
                 "max_num_sem", 
                 "start_age_advent",
                 "max_dist_advent",
                 "emission_rate_advent", 
                 "dmax_prop_advent", 
                 "max_num_advent",
                 "dmin",
                 "dmax",
                 "slope_elongation_diam", 
                 "type_tropism", 
                 "intensity_tropism", 
                 "coef_growth_duration", 
                 "inter_prim_distance", 
                 "proba_emergence_dmax", 
                 "proba_emergence_dmin", 
                 "ratio_diameters", 
                 "coef_variation_diam", 
                 "tissue_density", 
                 "growth_duration_coef", 
                 "life_expectancy_coef", 
                 "radial_growth_coef"
  )
  rs <- tibble(name = var_names, value = read_csv(path, col_names = F)$X1)%>% 
    spread(name, value)
  
  return(rs)
}

# Write the parameter file
write_archisimple <- function(rs, path){
  
  # names of the parameters
  var_names <- c("sim_length", 
                 "emission_rate_sem", 
                 "dmax_prop_sem", 
                 "max_num_sem", 
                 "start_age_advent",
                 "max_dist_advent",
                 "emission_rate_advent", 
                 "dmax_prop_advent", 
                 "max_num_advent",
                 "dmin",
                 "dmax",
                 "slope_elongation_diam", 
                 "type_tropism", 
                 "intensity_tropism", 
                 "coef_growth_duration", 
                 "inter_prim_distance", 
                 "proba_emergence_dmax", 
                 "proba_emergence_dmin", 
                 "ratio_diameters", 
                 "coef_variation_diam", 
                 "tissue_density", 
                 "growth_duration_coef", 
                 "life_expectancy_coef", 
                 "radial_growth_coef"
  )
  rs <- rs %>% 
    gather("name", "value")
  
  rs <- left_join(data.frame(name = var_names),    # Reorder data frame
            rs,
            by = "name")
  
  cat(rs$value, file=path, sep='\n') # Create the input file for Archisimple
}



