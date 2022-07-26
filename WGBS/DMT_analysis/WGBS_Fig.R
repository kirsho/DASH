#!/usr/bin/env Rscript

# WGBS_Fig.R v1.1 2022-07-13 with config file 
# Only for merged for now

# Take argument from sbatch
opts = commandArgs(trailingOnly=TRUE)

# recode in R
RDATA <- opts[1]  ## RData 
CONDI <- opts[2]  ## Conditions from the RData, ok when working from bams, must be rename for Figures
OUTDIR <- opts[3]
DESIGN <- opts[4]


CONDIORI <- CONDI
# Create outputfile in working directory
## Create path
REDIR2 <- gsub("script", OUTDIR, getwd())  
#REDIR3 <- REDIR2 # only for working from FigRData
## Initialise output file
sink(file = paste0(REDIR2, CONDI, "_Figures_output.txt"))

"#################################################################"
"Annotation and figures"
"https://www.bioconductor.org/packages/devel/bioc/vignettes/annotatr/inst/doc/annotatr-vignette.html"
"Anntotation with genomic features"
"Works with any genomic features, included or custom"
"CGI annotation"
"Genomic annotation"
"Use TxDb.* and org.*.eg.db packages"
"#################################################################"
cat("\n")
"#################################################################"
"Date"
start.time <- Sys.time()
start.time
"User : Olivier Kirsh"
"Mail : olivier.kirsh@u-paris.fr"
"Script : WGBS_Fig.R"
"opts -> arguments passed through sbatch"
opts
cat("\n")

"4 arguments"
"RData image"
RDATA
"Conditions"
CONDI
"Outdir"
OUTDIR
"Design"
DESIGN 

cat("\n")
### Load RData image
"#################################################################"
"Load RData image"
load(RDATA) 
"#################################################################"
cat("\n")

"#################################################################"
"Load personnal settings"
source("Methylator_Config.txt")
"#################################################################"
cat("\n")

"#################################################################"
"Check arguments"
"#################################################################"
cat("\n")

if(exists("REDIR")){
    print("REDIR")
    print(REDIR)
    print("is replaced by")
    print(REDIR2) # only for working from FigRData put back REDIR2
}

# set or erase REDIR path
REDIR <- REDIR2 # only for working from FigRData
cat("\n")
print("Path to working directory")
REDIR
cat("\n")
cat("\n")
print("Methylator_Config user arguments")
myarg

cat("\n")
"#################################################################"
"load packages"
"#################################################################"
library(xlsx)
library(data.table)
#library(mixtools)
#library(tzdb)
#library(readr)

library(magrittr)
library(dplyr)
library(ggplot2)
library(gridExtra)

library(GenomicRanges)
library(methylKit)
library(annotatr)
library(ggalluvial)
#library(edmr)
#library(genomation)



cat("\n")
"#################################################################"
"recover Sample ID"
"#################################################################"

if(exists("fu.tiles")){
    sample.ids <- myDiff.tiles@sample.ids
}else{
    sample.ids <- myDiff.CpG@sample.ids
}

sample.ids
cat("\n")

"sample.ids as vector ANA"
ANA <- paste0(sample.ids[2],"vs", sample.ids[1])


###############################################################################################
# For 2i vs Serum
#sample.ids[1] <- "Serum"
#sample.ids[2] <- "2i"

# For KO vs Serum
sample.ids[1] <- "Serum"
sample.ids[2] <- paste0(sample.ids[2], " (KO)")

###############################################################################################

"sample.ids as vector ANAtitle"
ANAtitle <- paste0(sample.ids[2]," vs ", sample.ids[1])

cat("\n")
"# load prepared Annots"
"Organism"
ORG
"Custom tracks"
CUSTOM

if(CUSTOM != "yes" & ORG == "mm10"){
load("Annotatr_mm10.RData")
}else{
 print("load Annotatr_mm10_Customs.RData")   
 load("Annotatr_mm10_Customs.RData")    # change for custom sans s when redone
}

cat("\n")
cat("annotations")
annotations

cat("\n")
cat("shortcuts")
shortcuts


cat("\n")
"Preparing dataframes for figures"
##  `AllDiffTiles` object definition

if(exists("fu.tiles")){
    methylKit::getData(fu.tiles) %>% dplyr::mutate(Methyl.Ctrl = numCs1/coverage1) -> fu.tilepc
    fu.tilepc %>% dplyr::mutate(Methyl.Cond = numCs2/coverage2) -> fu.tilepc
    dplyr::inner_join(myDiff.tiles %>% data.frame, fu.tilepc) -> AllDiffTiles
}else{
    methylKit::getData(fu.meth) %>% dplyr::mutate(Methyl.Ctrl = numCs1/coverage1) -> fu.tilepc
    fu.tilepc %>% dplyr::mutate(Methyl.Cond = numCs2/coverage2) -> fu.tilepc
    dplyr::inner_join(myDiff.CpG %>% data.frame, fu.tilepc) -> AllDiffTiles
}


cat("\n")
##  `AllDiffTiles` Add  Sample DNA methylation segmentation
### On Ctrl  `Methstatus_Ct`
AllDiffTiles <- AllDiffTiles %>% 
                      dplyr::mutate( Methstatus_Ct = ifelse(AllDiffTiles$Methyl.Ctrl <= 1/4, "Low",
                                                     ifelse(AllDiffTiles$Methyl.Ctrl >= 3/4, "High", "Mid")
                                                            )
                                    )

AllDiffTiles$Methstatus_Ct <- factor(AllDiffTiles$Methstatus_Ct, levels= c("Low", "Mid", "High"))

cat("\n")
### On Condition `Methstatus_Cd`
AllDiffTiles <- AllDiffTiles %>% 
                      dplyr::mutate( Methstatus_Cd = ifelse(AllDiffTiles$Methyl.Cond <= 1/4, "Low",
                                                     ifelse(AllDiffTiles$Methyl.Cond >= 3/4, "High", "Mid")
                                                            )
                                    )

AllDiffTiles$Methstatus_Cd <- factor(AllDiffTiles$Methstatus_Cd, levels= c("Low", "Mid", "High"))

cat("\n")
##  `AllDiffTiles` Add  Differential DNA methylation status `Diff_expr`
cat("\n")
DIF
cat("\n")
QV
cat("\n")
AllDiffTiles <- AllDiffTiles %>% dplyr::mutate(Diff_expr = ifelse(abs(meth.diff) >= DIF & qvalue <= QV,
                                           "Significant",
                                           "non-significant") 
                                        )

AllDiffTiles$Diff_expr <- factor(AllDiffTiles$Diff_expr, levels =c("Significant", "non-significant" ) )



##  `AllDiffTiles` Add  Differential DNA methylation status `DM_status`
AllDiffTiles <- AllDiffTiles %>% 
                        dplyr::mutate(DM_status =   ifelse(AllDiffTiles$meth.diff <= -DIF & AllDiffTiles$qvalue <= QV, "Hypo",
                                                    ifelse(AllDiffTiles$meth.diff >= DIF & AllDiffTiles$qvalue <= QV, "Hyper", "None")
                                                            )
                                     )

AllDiffTiles$DM_status <- factor(AllDiffTiles$DM_status, levels= c("Hyper", "Hypo", "None"))

cat("\n")
##  Check `AllDiffTiles`
cat("\n")
AllDiffTiles  %>%  head
cat("\n")
AllDiffTiles  %>%  nrow



## Creation of `dfpc` data.frame, vectorized version of `AllDiffTiles` 

dfpc <- data.frame( chr = rep( AllDiffTiles$chr, 2),
                    start = rep( AllDiffTiles$start, 2),
                    end = rep( AllDiffTiles$end, 2),
                    strand = rep( AllDiffTiles$strand, 2),
                    Meth_perc = c( AllDiffTiles$Methyl.Ctrl, AllDiffTiles$Methyl.Cond ),
                    Coverage = c( AllDiffTiles$coverage1, AllDiffTiles$coverage2 ),
                    ID = factor( rep( c( sample.ids[1], sample.ids[2] ),
                                 each = nrow( AllDiffTiles ) ),
                                 levels = c( sample.ids[1], sample.ids[2] ) ),
                    Diff_expr = c( AllDiffTiles$Diff_expr, AllDiffTiles$Diff_expr ),
                    DM_status = c( AllDiffTiles$DM_status, AllDiffTiles$DM_status ),
                    Methstatus = c( AllDiffTiles$Methstatus_Ct, AllDiffTiles$Methstatus_Cd ),
                    Methstatus_Ct = c( AllDiffTiles$Methstatus_Ct, AllDiffTiles$Methstatus_Ct ),
                    Methstatus_Cd = c( AllDiffTiles$Methstatus_Cd, AllDiffTiles$Methstatus_Cd )
                   )
cat("\n")
dfpc %>% head  
cat("\n")    
dfpc %>% colnames

cat("\n")
# region to annotate as GRanges
## `AllDiffTiles` as Granges
(AllDiffTiles_regions <- makeGRangesFromDataFrame( AllDiffTiles , 
                                        keep.extra.columns=TRUE,
                                        seqinfo=Seqinfo(genome=ORG),
                                        ignore.strand=F)
)

cat("\n")
## `dfpc` as Granges
(dfpc_regions <- makeGRangesFromDataFrame( dfpc , 
                                        keep.extra.columns=TRUE,
                                        seqinfo=Seqinfo(genome=ORG),
                                        ignore.strand=F)
)

cat("\n")
# Check annotations
(annotations)

cat("\n")
## `AllDiffTiles` annotated
(AllDiffTiles_annotated <- annotate_regions( regions = AllDiffTiles_regions,
                                            annotations = annotations,
                                            ignore.strand = TRUE,
                                            quiet = FALSE)
 )

cat("\n")
AllDiffTiles_annotated %>% data.frame() %>% head
cat("\n")
AllDiffTiles_annotated %>% data.frame() %>% colnames

cat("\n")
## `dfpc` annotated
(dfpc_annotated <- annotate_regions( regions = dfpc_regions,
                                            annotations = annotations,
                                            ignore.strand = TRUE,
                                            quiet = FALSE)
 )

cat("\n")
dfpc_annotated %>% data.frame() %>% head
cat("\n")
dfpc_annotated %>% data.frame() %>% colnames

cat("\n")
# Randomisation
## `AllDiffTiles_regions`
(random_AllDiffTiles_regions <- randomize_regions(regions = AllDiffTiles_regions,
                                                 allow.overlaps = TRUE,
                                                 per.chromosome = TRUE)
 )

cat("\n")
## `dfpc_regions`
 (random_dfpc_regions <- randomize_regions(regions = dfpc_regions,
                                                 allow.overlaps = TRUE,
                                                 per.chromosome = TRUE)
 
 )

cat("\n")
# Annotation of Randomised regions
## Annotate randomize region All
 (random_AllDiffTiles_annotated <- annotate_regions(regions = random_AllDiffTiles_regions,
                                                          annotations = annotations,
                                                          ignore.strand = TRUE,
                                                          quiet = TRUE) 
)

cat("\n")
## Annotate randomize region dfpc
(random_dfpc_annotated <- annotate_regions(regions = random_dfpc_regions,
                                                  annotations = annotations,
                                                  ignore.strand = TRUE,
                                                  quiet = TRUE)
)



######################################################################################
# Filter AlldiffTiles and make Signif GR
SigDiffTiles_region <- makeGRangesFromDataFrame(as.data.frame(AllDiffTiles) %>%
                                                    filter(., Diff_expr == "Significant"), 
                                                    keep.extra.columns=TRUE,
                                                    seqinfo=Seqinfo(genome="mm10"),
                                                    ignore.strand=F)
# Annot signif GR
SigDiffTiles_annotated <- annotate_regions(regions = SigDiffTiles_region,
                                                          annotations = annotations,
                                                          ignore.strand = TRUE,
                                                          quiet = TRUE)

# randomize Signif GR
random_SigDiffTiles_region <- randomize_regions(regions = SigDiffTiles_region,
                                                    allow.overlaps = TRUE,
                                                    per.chromosome = TRUE)

# Annot signif GR
random_SigDiffTiles_annotated <- annotate_regions(regions = random_SigDiffTiles_region,
                                                          annotations = annotations,
                                                          ignore.strand = TRUE,
                                                          quiet = TRUE)



sig <- as.data.frame(AllDiffTiles) %>% filter(., Diff_expr == "Significant")


dfpc_Sig <- data.frame( chr = rep( sig$chr, 2),
                    start = rep( sig$start, 2),
                    end = rep( sig$end, 2),
                    strand = rep( sig$strand, 2),
                    Meth_perc = c( sig$Methyl.Ctrl, sig$Methyl.Cond ),
                    Coverage = c( sig$coverage1, sig$coverage2 ),
                    ID = factor( rep( sample.ids , each = nrow( sig ) ), levels = sample.ids ),
                    Diff_expr = c( sig$Diff_expr, sig$Diff_expr ),
                    DM_status = c( sig$DM_status, sig$DM_status ),
                    Methstatus = c( sig$Methstatus_Ct, sig$Methstatus_Cd ),
                    Methstatus_Ct = c( sig$Methstatus_Ct, sig$Methstatus_Ct ),
                    Methstatus_Cd = c( sig$Methstatus_Cd, sig$Methstatus_Cd )
                   )
  

dfpc_Sig %>% head 

dfpc_Sig %>% colnames




# Filter AlldiffTiles and make Signif GR
dfpc_Sig_region <- makeGRangesFromDataFrame(as.data.frame(dfpc_Sig) %>%
                                                    filter(., Diff_expr == "Significant"), 
                                                    keep.extra.columns=TRUE,
                                                    seqinfo=Seqinfo(genome=ORG),
                                                    ignore.strand=F)
# Annot signif GR
dfpc_Sig_annotated <- annotate_regions(regions = dfpc_Sig_region,
                                                          annotations = annotations,
                                                          ignore.strand = TRUE,
                                                          quiet = TRUE)



## Creation of `methst` data.frame  vectorized version of `AllDiffTiles$Methstatus_Cd` &  `AllDiffTiles$Methstatus_Ct` ~ `AllDiffTiles$ID`

cat("\n")
"methst object"

methst <- data.frame( Methstatus = c( AllDiffTiles$Methstatus_Ct, AllDiffTiles$Methstatus_Cd ),
                      ID = factor( rep( c( sample.ids[1], sample.ids[2] ), each = nrow( AllDiffTiles ) ),
                                   levels=c( sample.ids[1], sample.ids[2] ) )
                      )

cat("\n")
"methst object" %>% print
methst %>% head  
methst %>% colnames


if(CUSTOMPLOTS != "yes" ){
 source("basic_plot.R")
}else{
 print("je source custom_plot.R")  
 source("custom_plot.R")
}

print("custom_plot.R runs well") 


print( " report SessionInfo" )
sessionInfo()

cat("\n")
"#################################################################"
cat("\n")

"An RData image is save in ../RData folder"
if(SAVERDATA == "yes" ){
 save.image(file =paste0("../RData/", "Figure_", CONDIORI,"_.RData") ) 
}else{
 print("it is done")
}


# Close output file
sink()