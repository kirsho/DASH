Olivier Kirsh (olivier.kirsh@u-paris.fr)  
2022-07-27



methylator.simg 
====  

These 3 files allow you to recreate a singularity image (methylator.simg) used in WGBS secondary analysis (DMC, DMT, eDMR, annotation and figures)
- methylator.yaml can be used to create a conda env (list of tools and R packages)    
- edmr_install.R is called through Rscript at the singularity build step (Installed via github)  
- Singularity : Singularity definition file  
