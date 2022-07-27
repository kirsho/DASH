#!/usr/bin/env Rscript

library(devtools)

install_github("ShengLi/edmr", 
                lib = "/opt/conda/envs/methylator/lib/R/library", 
                upgrade = "never", 
                force = 'yes')




