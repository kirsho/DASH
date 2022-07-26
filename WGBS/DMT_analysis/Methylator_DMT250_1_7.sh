#!/bin/bash
################################ Slurm options #################################

### Job name
#SBATCH --job-name=MergedDMTFig250-1-7_220722                    # ****Change****

### Output
##SBATCH --output=          # both STDOUT and STDERR 
#SBATCH -o ../outerr/MergedDMTFig250-1-7_220722.%A.%a.out         # STDOUT file with the Node name, the Job ID ****Change****
#SBATCH -e ../outerr/MergedDMTFig250-1-7_220722.%A.%a.err         # STDERR file with the Node name and the Job ID ****Change****

### Limit run time "days-hours:minutes:seconds"
##SBATCH --time=24:00:00  (max with fast partition = 24h, but you can put less)

### Requirements
#SBATCH --partition=ipop-up		                        # set to project name on ifbcore cluster 
#SBATCH --array=1                                     # define number of task (how many RData?) ****Change****

### Advanced options
#SBATCH --cpus-per-task=1                               # define number of cpu per task (leave to one if you don't know)
#SBATCH --mem=116G                                      # define amount of ram per task (Book enough but not too much (max is 256)

### More advanced options
##SBATCH --ntasks-per-node=1
##SBATCH --mem-per-cpu=5GB
##SBATCH --nodes=1

################################################################################
## Version 2022-06-28 V1.0 with config file
################################################################################

echo '########################################'
echo 'Date:' $(date --iso-8601=seconds)
echo 'User:' $USER
echo 'Host:' $HOSTNAME
echo 'Job Name:' $SLURM_JOB_NAME
echo 'Job Id:' $SLURM_JOB_ID
echo 'Directory:' $(pwd)
echo '########################################'
echo 'Methylator.sh'
echo 'Projet : DROM'                                                                # ****Change****
echo 'Analysis : Annotation & Figures WGBS DROM Serum, 2i & KOs merged OK'          # ****Change****
echo 'Mail : olivier.kirsh@u-paris.fr'                                              # ****Change****

start0=`date +%s`

################################################################################
### Choose analysis type
################################################################################

DMC=no            # explo or analysis or no
DMT=no            # analysis or no
EDMR=no           # yes or no
BEDG=no           # explo or diff
FIG=yes           # yes or no


################################################################################
### Choose experimental design  
################################################################################
# Examples
#### RData from count tables 
# methexport
# allmethexport
#
#### RData from bam
#DESIGN=merged
#DESIGN=replicates
#DESIGN=all

DESIGN=merged                                                                       #dev todo chang

################################################################################
### Defines INPUT Files
### Don't forget to define/change  the slurm option --array=1 --array=1-7 --array=1-100%10
### **** Change and/or remove # before desired option ****
################################################################################
#
### Methods1 / Give the path the file to with your MethylKit MethylRaw object.
### For count tables
# *-methexport.RData
### For design with merged
# Merged*flt.RData
### For design with replicates
# Replicates*flt.RData
### For analysis with multiple conditions
# All*flt.RData
### For figures
#  MethylKitDMT_Merged*vsS_flt_250_1_7_delta10_q5pc.RData

#RDATA=Merged*flt.RData     #dev todo change for cat .....$2 ...

RDATA=MethylKitDMT_MergedBEND3vsS_flt_250_1_7_delta10_q5pc.RData

################################################################################################################################################################
### don't touch anything bellow
################################################################################################################################################################
#
#             uu$$$$$$$$$$$uu
#           uu$$$$$$$$$$$$$$$$$uu
#         u$$$$$$$$$$$$$$$$$$$$$u
#        u$$$$$$$$$$$$$$$$$$$$$$$u
#       u$$$$$$$$$$$$$$$$$$$$$$$$$u
#       u$$$$$$$$$$$$$$$$$$$$$$$$$u
#       u$$$$$$"   "$$$"   "$$$$$$u
#       "$$$$"      u$u       $$$$"
#        $$$u       u$u       u$$$
#        $$$u      u$$$u      u$$$
#         "$$$$uu$$$   $$$uu$$$$"
#          "$$$$$$$"   "$$$$$$$"
#            u$$$$$$$u$$$$$$$u
#             u$"$"$"$"$"$"$u
#  uuu        $$u$ $ $ $ $u$$       uuu
# u$$$$        $$$$$u$u$u$$$       u$$$$
#  $$$$$uu      "$$$$$$$$$"     uu$$$$$$
#u$$$$$$$$$$$uu    """""    uuuu$$$$$$$$$$
#$$$$"""$$$$$$$$$$uuu   uu$$$$$$$$$"""$$$"
# """      ""$$$$$$$$$$$uu ""$"""
#           uuuu ""$$$$$$$$$$uuu
#  u$$$uuu$$$$$$$$$uu ""$$$$$$$$$$$uuu$$$
#  $$$$$$$$$$""""           ""$$$$$$$$$$$"
#   "$$$$$"                      ""$$$$""
#     $$$"                         $$$$"
#
#
################################################################################################################################################################

################################################################################
### Set output and results directory names 
################################################################################

INPUT=$(ls -A ${PWD%%script}RData/${RDATA} | sed -n ${SLURM_ARRAY_TASK_ID}p)  


################################################################################
### Extract Condition names
################################################################################

COND=$(basename ${INPUT} .RData) 


################################################################################
### Set output and results directory names 
################################################################################

### Exploration
OUTDIRcex=DMR${DESIGN}/methylkitDMC_Explo/
OUTDIRtex=DMR${DESIGN}/methylkitDMT_Explo/
OUTDIRb=DMR${DESIGN}/Bedgraph_Explo/                # remove Bedgraphcov${MINCOV}/, add value un file

### Analyses
OUTDIRmc=DMR${DESIGN}/methylkitDMC/
OUTDIRmt=DMR${DESIGN}/methylkitDMT/
OUTDIRe=DMR${DESIGN}/eDMR/

### Figures
OUTDIRfig=DMR${DESIGN}/${COND}_Figures/

################################################################################
### Create output directories 
################################################################################

### Exploration
[[ -d ${PWD%%script}${OUTDIRcex} ]] ||  mkdir -p ${PWD%%script}${OUTDIRcex}
[[ -d ${PWD%%script}${OUTDIRtex} ]] ||  mkdir -p ${PWD%%script}${OUTDIRtex}
[[ -d ${PWD%%script}${OUTDIRb} ]]   ||  mkdir -p ${PWD%%script}${OUTDIRb}

### Analyses
[[ -d ${PWD%%script}${OUTDIRmc} ]] ||  mkdir -p ${PWD%%script}${OUTDIRmc}
[[ -d ${PWD%%script}${OUTDIRmt} ]] ||  mkdir -p ${PWD%%script}${OUTDIRmt}
[[ -d ${PWD%%script}${OUTDIRe}  ]] ||  mkdir -p ${PWD%%script}${OUTDIRe}

### Figures
[[ -d ${PWD%%script}${OUTDIRfig} ]] ||  mkdir -p ${PWD%%script}${OUTDIRfig}


### Create a sub folder in ./RData  *** requires to update PATH un R script, not recommanded ***
# mkdir -p ${PWD%%script}RData/${DESIGN}/


################################################################################
### Deactivate any module
################################################################################

module purge

################################################################################
### Echo Variables in out/err files
################################################################################

echo "Echo Variables in outerr files"
echo "Input"
echo ${INPUT}
echo "Conditions"
echo ${COND}
echo "Design"
echo ${DESIGN}

echo "Jobs"
if [[ "$DMC" == "explo" ]]; then
echo "make DMC explo"
echo ${DMC}
echo "Outdir"
echo ${OUTDIRcex}
fi

if [[ "$DMC" == "analysis" ]]; then
echo "make DMC analysis"
echo ${DMC}
echo "Outdir"
echo ${OUTDIRmc}
fi

if [[ "$DMT" == "explo" ]]; then
echo "make DMT explo"
echo ${DMT}
echo "Outdir"
echo ${OUTDIRtex}
fi

if [[ "$DMT" == "analysis" ]]; then
echo "make DMT analysis"
echo ${DMT}
echo "Outdir"
echo ${OUTDIRmt}
fi

if [[ "$EDMR" == "yes" ]]; then
echo "make eDMR analysis"
echo ${EDMR}
echo "Outdir"
echo ${OUTDIRe}
fi

if [[ "$BEDG" == "explo" ]]; then
echo "make Bedgraph explo"
echo ${EDMR}
echo "Outdir"
echo ${OUTDIRb}
fi

if [[ "$BEDG" == "diff" ]]; then
echo "make Bedgraph diff"
echo ${EDMR}
echo "Outdir"
echo ${OUTDIRb}
fi

if [[ "$FIG" == "yes" ]]; then
echo "make figures"
echo ${FIG}
echo "Outdir"
echo ${OUTDIRfig}
fi



################################################################################
### Run R scripts 
# -B $PWD:/script
################################################################################

### DMC ################################################################################
if [[ "$DMC" == "explo" ]]; then
srun singularity exec ../methylator.simg Rscript --vanilla MKitDMC_Exploration.R ${INPUT} ${COND} ${OUTDIRcex} ${DESIGN}
fi 

if [[ "$DMC" == "analysis" ]]; then
srun singularity exec ../methylator.simg Rscript --vanilla MKitDMC_Analysis.R ${INPUT} ${COND} ${OUTDIRmc} ${DESIGN} #${MINCOV} ${DIFC} ${QVC}
fi
################################################################################


### DMT ################################################################################
if [[ "$DMT" == "explo" ]]; then
srun singularity exec ../methylator.simg Rscript --vanilla MKitDMT_Exploration.R ${INPUT} ${COND} ${OUTDIRtex} ${DESIGN}
fi 

if [[ "$DMT" == "analysis" ]]; then
srun singularity exec ../methylator.simg Rscript --vanilla MKitDMT_Analysis.R ${INPUT} ${COND} ${OUTDIRmt} ${DESIGN} # ${MINCOV} ${DIFC} ${QVC} ${TS} ${RTS} ${DIF} ${QV}
fi
################################################################################


### EDMR ################################################################################
if [[ "$EDMR" == "yes" ]]; then
srun singularity exec ../methylator.simg Rscript --vanilla MKit_EDMR_Analysis.R ${INPUT} ${COND} ${OUTDIRe} ${DESIGN} # ${MINCOV} ${DIFC} ${QVC} ${CPG} ${DMC} ${DIF} ${QV}
fi
################################################################################


### BEDGRAPH ################################################################################
if [[ "$BEDG" == "explo" ]]; then
srun singularity exec ../methylator.simg Rscript --vanilla MKit_Bedgraph.R ${INPUT} ${COND} ${OUTDIRb} ${DESIGN} ${MINCOV}
fi

if [[ "$BEDG" == "diff" ]]; then
srun singularity exec ../methylator.simg Rscript --vanilla MKit_BedgraphDiff.R ${INPUT} ${COND} ${OUTDIRb} ${DESIGN} ${MINCOV}
fi
################################################################################

### FIGURES ################################################################################
if [[ "$FIG" == "yes" ]]; then
srun singularity exec ../methylator.simg Rscript --vanilla WGBS_Fig.R ${INPUT} ${COND} ${OUTDIRfig} ${DESIGN} 
fi
################################################################################
