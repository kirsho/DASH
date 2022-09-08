Olivier Kirsh (olivier.kirsh@u-paris.fr)  
2022-07-27



methylator.simg 
====  

These 3 files allow you to recreate a singularity image (methylator.simg) used in WGBS secondary analysis (DMC, DMT, eDMR, annotation and figures)
   
- edmr_install.R is called through Rscript at the singularity build step (Installed via github).    
- Singularity : Singularity definition file to build a complete methylator Singularity image.    

- methylator.yaml can be used to create a conda env (list of tools and R packages) without eDMR.  It can be used easily with (yml2sing.sh)[https://github.com/kirsho/yml2sing] script  
