Olivier Kirsh (olivier.kirsh@u-paris.fr)  
2022-07-27



methylator.simg 
====  

These 2 following files allow you to recreate a singularity image (methylator.simg) used in WGBS secondary analysis (DMC, DMT, eDMR, annotation and figures)
   
- `edmr_install.R` is called through Rscript at the singularity build step (eDMR is distributed via github and can not be installed with conda).    
- `Singularity` : Singularity definition file to build a complete methylator Singularity image.    


Former Stuff:
`methylator.yaml` can be used to create a conda env (list of tools and R packages) without eDMR.  It can be used easily with (yml2sing.sh)[https://github.com/kirsho/yml2sing] script  
