
# Methylator package
Methylator is a bunch of R scripts designed for Methyl-Seq secondary analysis on HPC infrasctructure like ipop-up or ifb with **slurm** scheduler.  
Primary analysis must be done with **WGBS-workflow** which starts with `.fastq` files and outputs `.bam` files, `.bedgraph` files for nucleotides methylation levels and coverage and global methylation level. 

## Methylator can do :  
- DMC : Differentially methylated CpG analysis.  
- DMT : Differentially methylated Tiles analysis.  
- eDMR : Differentially methylated Region analysis.  
- BEDG : Generate `.bedgraph` files (% methylation, coverage, differential % methylation and q-value)  
- FIG : Annotation and figures.  

Methylator runs with on a Singularity/Apptainer container to respect the FAIR principes. You analysis is controlled by `Methylator_xxx.sh` to define and choose your analysis. A `Methylator_Config.txt` to pass specific parameters to your analysis. <!-- I wonder if some of the information in `Methylator_xxx.sh` should be given in `Methylator_Config.txt` -->

- `Methylator_xxx.sh`   
- `Methylator_Config.txt`   
- `MKitxxx_Analysis.R`  
  - `MKitDMC_Analysis.R`  
  - `MKitDMT_Analysis.R`  
  - `MKit_EDMR_Analysis.R`  
- `MKit_Bedgraph.R`  
- `MKit_BedgraphDiff.R`   
- `methylator.yaml`  
- `WGBS_Fig.R`  
- `custom_plot.R`  
- `Annotatr_mm10.R`  
- `Annotatr_mm10.RData`  
- `Annotatr_mm10_Custom_V2.R`  
- `Annotatr_mm10_Customs.RData`  

<!-- figures and annotatr stuff here C:\Users\olivi\Nextcloud\Bioinfo_Analyses\WGBSDROM_script\script\figures_V2  -->
<!-- annotation stuff here too C:\Users\olivi\Nextcloud\Bioinfo_Scripts_Tools\Annotations\Annot_with_mm10 -->


MKit_EDMR_Analysis.R




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
