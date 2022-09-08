
# Methylator package
Methylator is a bunch of R scripts designed for Methyl-Seq secondary analysis, annotations and figure on HPC infrasctructure like ipop-up or ifb with **slurm** scheduler.  
Primary analysis must be done with **WGBS-workflow** which starts with `.fastq` files and outputs `.bam` files, `.bedgraph` files for nucleotides methylation levels and coverage and global methylation level and an RData with BS in a (methylKit)[https://github.com/al2na/methylKit] object. 

## Methylator can do :  
- DMC : Differentially methylated CpG analysis.  
- DMT : Differentially methylated Tiles analysis.  
- eDMR : Differentially methylated Region analysis.  
- BEDG : Generate `.bedgraph` files (% methylation, coverage, differential % methylation and q-value)  
- FIG : Annotations and figures.  

Methylator runs with on a Singularity/Apptainer container to respect the FAIR principes. You analysis is controlled by `Methylator_xxx.sh` to choose your analysis . A `Methylator_Config.txt` must be edited to pass specific parameters to your analysis. A `Annotatr_mm10_Customs.RData` can also be provided for user specific annotation. A `custom_plot.R` can also be given for user-specific plots.  
<!-- I wonder if some of the information in `Methylator_xxx.sh` should be given in `Methylator_Config.txt` -->

### Methylator scripts & files  
#### Files to edit  
- `Methylator_xxx.sh`   
- `Methylator_Config.txt`  
- `custom_plot.R`   
- `Annotatr_mm10_Custom_V2.R` 
- 
#### Scripts (do not edit!)
 
- `MKitDMC_Analysis.R`  
- `MKitDMT_Analysis.R`  
- `MKit_EDMR_Analysis.R`  
- `MKit_Bedgraph.R`  
- `MKit_BedgraphDiff.R`   
- `WGBS_Fig.R`  
- `plot.R` (this file doesn't exist yet)  
- `Annotatr_mm10.R` 
- `Annotatr_mm10.RData`  
- `Annotatr_mm10_Customs.RData`   

#### Singularity image
- `methylator.simg'

<!-- figures and annotatr stuff here C:\Users\olivi\Nextcloud\Bioinfo_Analyses\WGBSDROM_script\script\figures_V2  -->
<!-- annotation stuff here too C:\Users\olivi\Nextcloud\Bioinfo_Scripts_Tools\Annotations\Annot_with_mm10 -->


### Methylator Singularity container
All Apps, R and R packages are containerized in a Singularity image. This allows a strict version control and reproducibility of the analyses.  Please read the Singularity documentation or my (digest)[https://github.com/kirsho/Singularity/blob/master/Intro2Singularity.md] to learn more about Singularity.  
The `Singularity_Methylator` image can be recreated with this 2 files stored (here)[https://github.com/kirsho/DASH/tree/main/IMAGES/methylator] :  
- `Singularity`      
- `edmr_install.R`  






- Methylator_xxx.sh
sbatch file for WGBS secondary analysis.

You can choose which analysis you want to perfo


```{R}
################################################################################
### Choose analysis type
################################################################################

DMC=no            # explo or analysis or no
DMT=no            # analysis or no
EDMR=no           # yes or no
BEDG=no           # explo or diff
FIG=yes           # yes or no
```
