
## IO functions to read and write archisimple parameter files
library(tidyverse)
library(xml2)
### READ function

path <- "~/Dropbox/science/projects/archisimple/archisimple/src/archisimple93/parameter2.xml"
options(scipen = 999) 

## READ AN ARCHISIMPLE XML FILE AND TRANFORM INTO A TABLE
read_archisimple_xml <- function(path){
  
  xml <- read_xml(path)
  dat <- tibble(name = rep("", length(xml_children(xml))), value = rep(0, length(xml_children(xml))))
  
  for(i in c(1:length(xml_children(xml)))){
    
    dat$name[i] <-xml_child(xml, search = i) %>% xml_name()
    
    if(dat$name[i] != "exportName"){
      dat$value[i] <- xml_child(xml, search = dat$name[i]) %>% 
        xml_text() %>% 
        as.numeric()
    }else{
      dat$value[i] <- xml_child(xml, search = dat$name[i]) %>% 
        xml_text()
    }
  }
  
  dat <- dat %>% 
    spread(name, value)
  return(dat)
}


# Write the parameter file
write_archisimple_XML <- function(rs, name = "generic", species = "none", path){
  
  xml <- paste0("<Plant name = '",name,"' species = '",species,"'>")
  
  rs <- rs %>% gather("name", "value")
  
  for(i in c(1:nrow(rs))){
    
    xml <- paste0(xml, "<",rs$name[i],">",rs$value[i],"</",rs$name[i],">")  
  }
  xml <- paste0(xml, "</Plant>")
  xml2 <- read_xml(xml)
  
  write_xml(xml2, file=path) # Create the input file for Archisimple
}



