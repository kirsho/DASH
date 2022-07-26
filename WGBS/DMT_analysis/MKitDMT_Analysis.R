#!/usr/bin/env Rscript

# MKit_DMT_Analysis.R v3.0.9 2022-07-21 ("." in dec .bed)
# Compute DiffMeth on tiles
# http://compgenomr.github.io/book/bsseq.html
# https://bioconductor.org/packages/release/bioc/vignettes/methylKit/inst/doc/methylKit.html
# https://github.com/al2na/methylKit
# more arguments passed through methylator
# Developped & Maintained by Olivier Kirsh <olivier.kirsh@u-paris.fr> 
# perform MethylKit DMTiles analysis 
# launch as Rscript with Methylator.sh srun Rscript --vanilla MKitDMT_Analysis.R ${INPUT} ${COND} ${OUTDIR} ${DESIGN} 
# requires a RData environment with a methylKit object ex:  ConditionvsControl.RData


### Arguments list defined in Methylator_Config.txt

# Take argument from sbatch
opts = commandArgs(trailingOnly=TRUE)

# recode in R
RDATA <- opts[1]  ## RData 
CONDI <- opts[2]  ## Conditions from the RData, ok when working from bams, must be rename for Figures
OUTDIR <- opts[3]
DESIGN <- opts[4]

# Create outputfile in working directory
## Create path
REDIR <- gsub("script", OUTDIR, getwd())  
## Initialise output file

"#################################################################"
"Load personnal settings"
source("./Methylator_Config.txt")
"#################################################################"
cat("\n")

### Save console output in a file in OUTDIR directory
sink(file = paste0(REDIR,"MethylKitDMT_", CONDI,"_", TS, "_", RTS, "_", COV, "-output.txt"))

"#################################################################"
"MethylKit_DMT_Analysis"
"adapted from https://github.com/al2na/methylKit"
"#################################################################"
"Date"
start.time <- Sys.time()
start.time
"User : Olivier Kirsh"
"Mail : olivier.kirsh@u-paris.fr"
"Script : MKit_DMT_Analysis.R"


cat("\n")
### Load RData image
"#################################################################"
"Load RData image"
load(RDATA) 
cat("\n")

"#################################################################"
"Check arguments"
"#################################################################"
cat("\n")

"opts passed by Methylator & config file"
opts
cat("\n")
"RData image"
RDATA
"Conditions"
CONDI
"Outdir"
OUTDIR
"Design"
DESIGN 
"Coverage"
COV
"SigDiffMeth"
DIF
"Q Value"
QV
"Tiles Size"
TS
"Relative Step Size"
RTS

### load packages
"#################################################################"
"load packages"

library(xlsx)
#library(mixtools)
#library(data.table)
#library(tzdb)
#library(readr)

library(magrittr)
library(dplyr)
library(ggplot2)
library(gridExtra)

library(GenomicRanges)
library(methylKit)
#library(TxDb.Mmusculus.UCSC.mm10.knownGene)
#library(annotatr)
#library(edmr)
#library(genomation)

cat("\n")
"#################################################################"
"recover Sample ID"
"#################################################################"
for(i in 1:length(myobj)){
    print("sample ID") 
    print(myobj[[i]]@sample.id)
    cat("\n")
}

cat("\n")
"myobj"
myobj
cat("\n")

"#################################################################"
"MethylKit Tiling window approach, Create Tiles"
"#################################################################"
tiles = tileMethylCounts(myobj,
                         win.size = TS,
                         step.size = TS * RTS,
                         cov.bases = COV )
cat("\n")
"tiles obj"
tiles
cat("\n")

"#################################################################"
"Unite tiles among samples"
"#################################################################"

## Filtering COVERAGE and remove top 0.1%
f.tiles=filterByCoverage(tiles,
                        lo.perc=NULL,
                        hi.count=NULL,
                        hi.perc=99.9)

cat("\n")
"f.tiles obj"
f.tiles
cat("\n")


"#################################################################"
"Methylation & coverage Statistics on raw data (tiles) & filtered data (f.tile) "
"#################################################################"

cat("\n")

for(i in 1:length(tiles)){
    print("tiles, raw data")
    print("sample ID")
    print(tiles[[i]]@sample.id)
    cat("\n")
    getMethylationStats(tiles[[i]],
                            plot=FALSE,both.strands=FALSE)
    cat("\n")
    getCoverageStats(tiles[[i]],
                            plot=FALSE,both.strands=FALSE)
}

for(i in 1:length(f.tiles)){
    print("f.tiles, filtered data")
    print("sample ID")
    print(f.tiles[[i]]@sample.id)
    cat("\n")
    getMethylationStats(f.tiles[[i]],
                            plot=FALSE,both.strands=FALSE)
    cat("\n")
    getCoverageStats(f.tiles[[i]],
                            plot=FALSE,both.strands=FALSE)
}


"#################################################################"
"Unite tiles among samples"
"#################################################################"
fu.tiles=unite(f.tiles,
              destrand=FALSE, # keep the bases covered in all samples TRUE for CpG only and base paire resolution
              chunk.size = 1e+06,
              mc.cores = 1)  # set to 1 if working locally on a Windows PC
"u.tiles obj"
fu.tiles
cat("\n")


cat("\n")
"Total Numbers of tiles in filtered & united tiles"
dim(fu.tiles)
cat("\n")


"#################################################################"
"Compute differential methylation Condition vs Control"
"#################################################################"
myDiff.tiles=calculateDiffMeth(fu.tiles,
                                covariates = NULL,
                                overdispersion = "none" ,      # c("none", "MN", "shrinkMN"),
                                adjust =  "SLIM" ,             # c("SLIM", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none", "qvalue"),
                                effect = "wmean",              # c("wmean", "mean", "predicted"),
#                               test = opts[],                 #  opts[] c("F", "Chisq", "fast.fisher", "midPval"),
                                mc.cores = 1 ,                 # opts[] set to 1 on windwows PC
                                slim = TRUE,
                                weighted.mean = TRUE,
                                chunk.size = 1e+06 ) 
cat("\n")
"myDiff.tiles object"
myDiff.tiles
cat("\n")
"str(myDiff.tiles)"
str(myDiff.tiles)
cat("\n")


"#################################################################"
"Filtering myDiff.tiles"
"#################################################################"
"my Tiles size"
TS
"relative Step size"
RTS
"Delta Diff %"
DIF
"QValue"
QV

cat("\n")
"#################################################################"
"Filtering of myDiff.tiles with Delta Diff = DIF et QValue = QV"
"myDiff.tiles.sig"
"myDiff.tiles.sighypo"
"myDiff.tiles.sighyper"

cat("\n")
myDiff.tiles.sig=getMethylDiff(myDiff.tiles,difference=DIF,qvalue=QV,type="all")
myDiff.tiles.sighypo=getMethylDiff(myDiff.tiles,difference=DIF,qvalue=QV,type="hypo")
myDiff.tiles.sighyper=getMethylDiff(myDiff.tiles,difference=DIF,qvalue=QV,type="hyper")

cat("\n")

"Filtering of myDiff.tiles with generic values |Delta Diff| > 25% and 3 qvalues 5% 1% 0.1%"
cat("\n")
"myDiff.tiles.sig005"
"myDiff.tiles.sig005hypo"
"myDiff.tiles.sig005hyper"
cat("\n")
"myDiff.tiles.sig001"
"myDiff.tiles.sig001hypo"
"myDiff.tiles.sig001hyper"
cat("\n")
"myDiff.tiles.sig0001"
"myDiff.tiles.sig0001hypo"
"myDiff.tiles.sig0001hyper"

# filter for difference=25,qvalue=0.05
myDiff.tiles.sig005=getMethylDiff(myDiff.tiles,difference=25,qvalue=0.05,type="all")
myDiff.tiles.sig005hypo=getMethylDiff(myDiff.tiles,difference=25,qvalue=0.05, type="hypo")
myDiff.tiles.sig005hyper=getMethylDiff(myDiff.tiles,difference=25,qvalue=0.05, type="hyper")

# filter for difference=25,qvalue=0.01
myDiff.tiles.sig001=getMethylDiff(myDiff.tiles,difference=25,qvalue=0.01,type="all")
myDiff.tiles.sig001hypo=getMethylDiff(myDiff.tiles,difference=25,qvalue=0.01, type="hypo")
myDiff.tiles.sig001hyper=getMethylDiff(myDiff.tiles,difference=25,qvalue=0.01, type="hyper")

# filter for difference=25,qvalue=0.001
myDiff.tiles.sig0001=getMethylDiff(myDiff.tiles,difference=25,qvalue=0.001,type="all")
myDiff.tiles.sig0001hypo=getMethylDiff(myDiff.tiles,difference=25,qvalue=0.001, type="hypo")
myDiff.tiles.sig0001hyper=getMethylDiff(myDiff.tiles,difference=25,qvalue=0.001, type="hypo")


# Create list for Statistics output
TilesListsig <- list(
    myDiff.tiles.sig,
    myDiff.tiles.sighypo,
    myDiff.tiles.sighyper
)

TilesListsig5pc <- list(
    myDiff.tiles.sig005,
    myDiff.tiles.sig005hypo,
    myDiff.tiles.sig005hyper
)

TilesListsig1pc <- list(    
    myDiff.tiles.sig001,
    myDiff.tiles.sig001hypo,
    myDiff.tiles.sig001hyper
)    

TilesListsig01pc <- list(    
    myDiff.tiles.sig0001,
    myDiff.tiles.sig0001hypo,
    myDiff.tiles.sig0001hyper
)

names(TilesListsig) <- c(
    "sig_",
    "sig_hypo",
    "sig_hyper"
)

names(TilesListsig5pc) <- c(
    "sig_005",
    "sig_005hypo",
    "sig_005hyper"
)

 names(TilesListsig1pc) <- c(   
    "sig_001",
    "sig_001hypo",
    "sig_001hyper"
)

names(TilesListsig01pc) <- c(    
    "sig_0001",
    "sig_0001hypo",
    "sig_0001hyper"
)



cat("\n")
"#################################################################"
"Some Statistics on Significant, Hyper and Hypo Tiles"
"#################################################################"
cat("\n")
"Total Numbers of tiles in filtered & united tiles"
dim(fu.tiles)
cat("\n")

"Total Numbers of tiles in DiffMeth all no filter"
dim(myDiff.tiles)
cat("\n")

"Statistics on myDiff.tiles.sig with chosen values"
cat("\n")
"   Numbers of signficant diff tiles Condition vs Control"
sapply(TilesListsig,nrow)

"   % of signficant diff tiles Condition vs Control"
round(100*sapply(TilesListsig,nrow)/nrow(myDiff.tiles),2)
cat("\n")


"Statistics myDiff.tiles with fixed values of Diff meth 25 % and 3 Q Values"
cat("\n")

"Statistics on significant (qvalue 5%) diff tiles "
"   Numbers of signficant diff tiles Condition vs Control"
sapply(TilesListsig5pc,nrow)

"   % of signficant diff tiles Condition vs Control"
round(100*sapply(TilesListsig5pc,nrow)/nrow(myDiff.tiles),2)
cat("\n")

"Statistics on significant (qvalue 1%) diff tiles "
"   Numbers of signficant diff tiles Condition vs Control"
sapply(TilesListsig1pc,nrow)

"   % of signficant diff tiles Condition vs Control"
round(100*sapply(TilesListsig1pc,nrow)/nrow(myDiff.tiles),2)
cat("\n")

"Statistics on significant (qvalue 0.1%), diff tiles "
"   Numbers of signficant diff tiles Condition vs Control"
sapply(TilesListsig01pc,nrow)

"   % of signficant diff tiles Condition vs Control"
round(100*sapply(TilesListsig01pc,nrow)/nrow(myDiff.tiles),2)
cat("\n")
cat("\n")


# "#################################################################"
# "Stats & PLots diffMethPerChromosomes"
# "#################################################################"
# cat("\n") 

# paste0("All DMT per chromosome,  diff = 0, QV=1")
# diffMethPerChr(myDiff.tiles, plot=FALSE, qvalue.cutoff=1, meth.cutoff=0)
# cat("\n")
# cat("\n")

# paste0("Significant DMT per chromosome,  diff ",DIF,"%, Qvalue ",QV)
# diffMethPerChr(myDiff.tiles, plot=FALSE, qvalue.cutoff=QV, meth.cutoff=DIF)
# cat("\n")
# cat("\n")

# pdf( paste0(REDIR, "MethylKitDMT_ChromPlot_MethDiff_",CONDI,".pdf"))
# diffMethPerChr(myDiff.tiles,plot=TRUE,qvalue.cutoff=1, meth.cutoff=0)
# diffMethPerChr(myDiff.tiles.sig,plot=TRUE,qvalue.cutoff=QV, meth.cutoff=DIF)
# dev.off()

# " See MethylKitDMT_ChromPlot_MethDiff_CondvsCtrl.pdf, 2 pages"
# cat("\n")  
# cat("\n")


"#################################################################"
"Write bed files for vizualisation"
"#################################################################"
cat("\n")
#Header : chr	start   end	strand	pvalue	qvalue	meth.diff
# Write bed myDiff.tiles
write.table(getData(myDiff.tiles)[, c(1:3,7,6)], 
            paste0(REDIR, "MethylKitDMT_myDiff.tiles_",CONDI,"_",TS,"_",COV,".bed"),
            sep='\t',
            dec=".",
            col.names=F,
            row.names=F,
            quote=F)
"Write bed file -> Done : See MethylKitDMT_myDiff.tiles_XXX.bed"            
cat("\n")

# Write bed XvsS_delta DIF QV.bed"
write.table(getData(myDiff.tiles.sig)[, c(1:3,7,6)], 
            paste0(REDIR, "MethylKitDMT_",CONDI,"_mydelta", DIF, "q",QV*100,"pc.bed"),
            sep='\t',
            dec=".",
            col.names=F,
            row.names=F,
            quote=F)
"Write bed file -> Done : See MethylKitDMT_CondvsCtrl_mydeltaXXqXpc.bed"            
cat("\n")

# Write bed XvsS_delta25q005.bed"
write.table(getData(myDiff.tiles.sig005)[, c(1:3,7,6)], 
            paste0(REDIR, "MethylKitDMT_",CONDI,"_delta25q5pc.bed"),
            sep='\t',
            dec=".",
            col.names=F,
            row.names=F,
            quote=F)

"Write bed file -> Done : See MethylKitDMT_CondvsCtrl_delta25q5pc.bed"            
cat("\n")

# Write bed XvsS_delta25q001.bed"
write.table(getData(myDiff.tiles.sig001)[, c(1:3,7,6)], 
            paste0(REDIR, "MethylKitDMT_",CONDI,"_delta25q1pc.bed"),
            sep='\t',
            dec=".",
            col.names=F,
            row.names=F,
            quote=F)

"Write bed file -> Done : See MethylKitDMT_CondvsCtrl_delta25q1pc.bed"             
cat("\n")

# Write bed XvsS_delta25q0001.bed"            
write.table(getData(myDiff.tiles.sig0001)[, c(1:3,7,6)], 
            paste0(REDIR, "MethylKitDMT_",CONDI,"_delta25q01pc.bed"),
            sep='\t',
            dec=".",
            col.names=F,
            row.names=F,
            quote=F)

"Write bed file -> Done : See MethylKitDMT_CondvsCtrl_delta25q1pc.bed"             
cat("\n")            

"#################################################################"
"Write txt & XLS files for exploration"
"#################################################################"
cat("\n")
# Write All tiles in a text file
write.table(data.frame(myDiff.tiles), 
            paste0(REDIR, "MethylKitDMT_all_",CONDI,".txt.xls"),
            sep='\t',
            dec=",",
            col.names=F,
            row.names=F,
            quote=F)
"Write All Tiles -> Done : See MethylKitDMT_all_CondvsCtrl.txt.xls (!!! it's a txt file, more than 65530 lines)"  
"Header : chr	start   end	strand	pvalue	qvalue	meth.diff"

cat("\n")

# Write Significant Tiles"
if (data.frame(myDiff.tiles.sig005) %>% nrow() < 65530 & 
    data.frame(myDiff.tiles.sig) %>% nrow() < 65530         ){

write.xlsx(data.frame(myDiff.tiles.sig),
            paste0(REDIR, "MethylKitDMT_signif_",CONDI,".xls"),
            sheetName = paste0("Tiles_myDiff_delta", DIF, "q", QV*100,"pc"),
            col.names=T,
            row.names=F,
            append = F) 

write.xlsx(data.frame(myDiff.tiles.sig005),
            paste0(REDIR, "MethylKitDMT_signif_",CONDI,".xls"),
            sheetName = "Tiles_delta25q005", 
            col.names=T,
            row.names=F,
            append = T)                       

# Write bed XvsS_delta25q001.txt"
write.xlsx(data.frame(myDiff.tiles.sig001),
            paste0(REDIR, "MethylKitDMT_signif_",CONDI,".xls"),
            sheetName = "Tiles_delta25q001", 
            col.names=T,
            row.names=F,
            append = T)    

# Write bed XvsS_delta25q0001.txt"            
write.xlsx(data.frame(myDiff.tiles.sig0001),
            paste0(REDIR, "MethylKitDMT_signif_",CONDI,".xls"),
            sheetName = "myTiles_delta25q0001", 
            col.names=T,
            row.names=F,
            append = T)             

print ("Write Significant Tiles -> Done : See MethylKitDMT_signif_CondvsCtrl.xls" )
cat("\n")

}else{

write.table(data.frame(myDiff.tiles.sig0001), 
            paste0(REDIR, "MethylKitDMT_signif_Diff25q01pc_",CONDI,".txt.xls"),
            sep='\t',
            dec=",",
            col.names=F,
            row.names=F,
            quote=F)

write.table(data.frame(myDiff.tiles.sig),
            paste0(REDIR, "MethylKitDMT_myDiff_delta", DIF, "q", QV*100,"pc_", CONDI,".txt.xls"),
            sep='\t',
            dec=",",
            col.names=F,
            row.names=F,
            quote=F)  

print("MethylKitDMT_signif_CondvsCtrl.xls is not created because more than 65536 lines"  )  
cat("\n")
print("")
print("Write Significant Tiles Diff25QV01pc and your Values -> Done ")
print(": See MethylKitDMT_signif_Diff25q01pc_CondvsCtrl.txt.xls" )
print(": See MethylKitDMT_myDiff_deltaXXqXXpc_CondvsCtrl.txt.xls" )
}





"#################################################################"
"Vizualisation "
"#################################################################"

# P Value distribution
pvtiles = ggplot(data=myDiff.tiles,
                    aes(x=pvalue)) +
                    xlim(0, 1) +
                    geom_vline(xintercept=c(0.05), col="red") +
                    geom_histogram(color="black", fill="white") +
                    theme_classic() +
                    ggtitle(paste0( "Diff_Tiles_P_Value_Distribution_",CONDI,"\n",
                                    "Good if not uniformally distributed"))

pdf( paste0(REDIR, "MethylKitDMT_Pvalue_Distribution_MethDiff_",CONDI,".pdf")) 
    pvtiles
dev.off()

cat("\n")
"Plot myDiff.tiles pValue Histogram  -> Done : See MethylKitDMT_Pvalue_Distribution_MethDiff_CondvsCtrl.pdf" 
cat("\n")


# Q Value distribution
qvtiles = ggplot(data=myDiff.tiles,
                    aes(x=qvalue)) +
                    xlim(0, 1) +
                    geom_vline(xintercept=c(0.05), col="red") +
                    geom_histogram(color="black", fill="white") +
                    theme_classic() +
                    ggtitle(paste0( "Diff_Tiles_Q_Value_Distribution_",CONDI))

pdf( paste0(REDIR, "MethylKitDMT_Qvalue_Distribution_MethDiff_",CONDI,".pdf"))
    qvtiles
dev.off()

cat("\n")
"Plot myDiff.tiles qValue Histogram  -> Done : See MethylKitDMT_Qvalue_Distribution_MethDiff_CondvsCtrl.pdf"
cat("\n")

# Volcano Plot 
#(myDiff.tiles)[sample(nrow(myDiff.tiles), round(nrow(myDiff.tiles)*0.01)), ] 
vctiles = ggplot(data=myDiff.tiles,
                    aes(x=meth.diff, y=-log10(qvalue))) +
                    xlim(-101, 101) +
                    geom_vline(xintercept=c(-25, 25), col="red") +
                    geom_hline(yintercept=-log10(c(0.05, 0.01, 0.001)), col="red") +
                    geom_point() +
                    theme_classic() +
                    ggtitle(paste0("VolcanoPLot_TilesMethDiff_",CONDI))

jpeg( paste0(REDIR, "MethylKitDMT_VolcanoPlot_MethDiff_",CONDI,".jpg"))
    vctiles
dev.off()

"Volcano Plot myDiff.tiles -> Done : See MethylKitDMT_VolcanoPlot_MethDiff_CondvsCtrl.jpg"

cat("\n")
cat("\n")
cat("\n")


"#################################################################"
"Annotation & intersection with genomic features"
"#################################################################"
cat("\n")
"#################################################################"
"Select filter for All chromosomes (allC) or autosomes (autoC)"
allC = c(paste0("chr", 1:21), "chrX","chrY")
autoC = c(paste0("chr", 1:21))

"allC"
CHROM <- allC
CHROM

cat("\n")
"#################################################################"
"Reference Gene intervals"
"  - read ../Ref/mm10refseqann.bed"
mm10refseqann <- read.delim("../Ref/mm10refseqann.bed", header=FALSE)
dim(mm10refseqann)
head(mm10refseqann)
cat("\n")

"  - Filter Chromosome set in mm10refseqann"
mm10refseqann <- mm10refseqann[ mm10refseqann$V1 %in% CHROM, ]
dim(mm10refseqann)
head(mm10refseqann)
cat("\n")

"  - Set as GRanges object : mm10Genes.gr"
colnames(mm10refseqann)=c("chrom","Start","End","strand","symbol","ensembl","entrez","refseq")
mm10Genes.gr <-as(mm10refseqann,"GRanges")
mm10Genes.gr
cat("\n")

"#################################################################"
"Reference promoters, filter for -1500 +1000 around TSS : prom.gr "
cat("\n")
prom.gr = GenomicRanges::promoters(mm10Genes.gr,upstream=1500, downstream=1000, use.names=F)
prom.gr
cat("\n")


"#################################################################"
"Reference CGI intervals "
"  - read ../Ref/mm10_CpG_Island.bed.txt"
mm10_CpG_Island <- read.delim("../Ref/mm10_CpG_Island.bed.txt")
mm10_CpG_Island <- mm10_CpG_Island[,-1]
dim(mm10_CpG_Island)
head(mm10_CpG_Island)
cat("\n")

"  - Filter Chromosome set in mm10_CpG_Island.bed"
mm10_CpG_Island <- mm10_CpG_Island[ mm10_CpG_Island$chrom %in% CHROM, ]
dim(mm10_CpG_Island)
head(mm10_CpG_Island)
cat("\n")

"  - Set as GRanges object : mm10CGI.gr"
mm10CGI.gr <-as(mm10_CpG_Island,"GRanges")
mm10CGI.gr
cat("\n")

"#################################################################"
"Reference CGI annotated with the closest promoters -1500 +1000 : mm10CGIprom.gr"
cat("\n")

clo <- nearest(mm10CGI.gr, prom.gr)
dist2clo <- distanceToNearest(mm10CGI.gr, prom.gr)
mm10_CpG_Island_annot_prom <- data.frame(mm10CGI.gr, 
                                        Dist2_= dist2clo@elementMetadata@listData$distance,
                                        symbol = prom.gr[clo,]$symbol,
                                        ensembl = prom.gr[clo,]$ensembl,
                                        entrez = prom.gr[clo,]$entrez)

mm10CGIprom.gr <- as(mm10_CpG_Island_annot_prom,"GRanges")
mm10CGIprom.gr
cat("\n")
cat("\n")

"#################################################################"
"myDiff Tiles to annotate"
cat("\n")
"  - myDiff.tiles, all the analysed Tiles"
myDiff.tiles.gr <- as(myDiff.tiles,"GRanges")
cat("\n")

"  - myDiff.tiles.sig, Significant Tiles"
myDiff.tiles.sig.gr <- as(myDiff.tiles.sig,"GRanges")
cat("\n")

"  - myDiff.tiles.sighyper, Significant hyper meth Tiles"
myDiff.tiles.sighyper.gr <- as(myDiff.tiles.sighyper,"GRanges")
cat("\n")

"  - myDiff.tiles.sighypo, Significant hypo meth  Tiles"
myDiff.tiles.sighypo.gr <- as(myDiff.tiles.sighypo,"GRanges")
cat("\n")

"Creation of myDiff.tiles.grl, a list of 4 Genomic ranges"
myDiff.tiles.grl <- GRangesList(myDiff.tiles.gr,
                    myDiff.tiles.sig.gr,
                    myDiff.tiles.sighyper.gr,
                    myDiff.tiles.sighypo.gr
                    )
cat("\n")
cat("\n")

##################################################################################################################################
# Set ExportDiffAnnot() function. Worg with GRanges objects
# 3 Arguments : 
#       difmeth =>  myDiff.tiles.grl or myDiff.tiles.gr
#       annot   =>  mm10CGIprom.gr or mm10Gene.gr or prom.gr
#       type    =>  "CGI" or "GeneBody" or "Prom"

ExportDiffAnnot <- function(difmeth, annot, type){                                                          # 3 args
                            overlap <- subsetByOverlaps(difmeth, annot)                                     # Subset Diffmeth by overlap with annot
                            dist2annot <- distanceToNearest(overlap, annot)                                 # Calculate distance between features
                            exportdiffannot <- data.frame(overlap, Overlap2_ = annot[dist2annot@to])        # Export annoted Diff meth
                            return(exportdiffannot)
}

# report 1st element of an overlap
# Change with findoverlap developped in MethByRanges
##################################################################################################################################


"#################################################################"
"Extraction & Annotation of myDiff Tiles that overlap with CGI"
"#################################################################"
difmeth = myDiff.tiles.grl
annot = mm10CGIprom.gr
type = "CGI"

difmeth.inter.CGI <- lapply(difmeth, ExportDiffAnnot, annot, type)

names(difmeth.inter.CGI) <- c("all.overlap",
                              "Sig.overlap",
                              "Sig.Hyper.overlap",
                              "Sig.Hypo.overlap") 
                              
cat("\n")

"Total numbers of CGI"
length(mm10CGI.gr)
cat("\n")

"Total numbers of Tiles"
length(myDiff.tiles.gr)
cat("\n")

"numbers of Tiles overlapping with a CGI"
sapply(difmeth.inter.CGI, nrow)
cat("\n")

"% of Tiles overlapping with a CGI"
round(100 * sapply(difmeth.inter.CGI, nrow)/nrow(myDiff.tiles) , 2)
cat("\n")


"Head of Tiles overlapping with a CGI (GR List of 4 Genomic Ranges)"
lapply(difmeth.inter.CGI, head)
cat("\n")

"Writting MethylKitDMT_all.overlap_CGI "
for(i in 1:length(difmeth.inter.CGI) ) {
  write.table(difmeth.inter.CGI[[i]],
              paste0(REDIR, "MethylKitDMT_",
                     names(difmeth.inter.CGI[i]),
                     "_CGI_myDiff_delta", 
                     DIF, "q", QV*100,"pc_", CONDI,".tsv"),
              sep='\t',
              dec=",",
              col.names=T,
              row.names=F,
              quote=F) 
}

cat("\n")
" Check the 4 outputs as ..."
"MethylKitDMT_all.overlap_CGI_myDiff_xxxxxxxxxx.tsv"
"MethylKitDMT_Sig.overlap_CGI_myDiff_xxxxxxxxxx.tsv"
"MethylKitDMT_HyperSig.overlap_CGI_myDiff_xxxxxxxxxx.tsv"
"MethylKitDMT_HypoSig.overlap_CGI_myDiff_xxxxxxxxxx.tsv"
cat("\n")
cat("\n")

"#################################################################"
"Extraction & Annotation of myDiff Tiles that overlap with Promoters"
"#################################################################"
difmeth = myDiff.tiles.grl
annot = prom.gr
type = "Prom"

difmeth.inter.PROM <- lapply(difmeth, ExportDiffAnnot, annot, type)

names(difmeth.inter.PROM) <- c("all.overlap",
                              "Sig.overlap",
                              "Sig.Hyper.overlap",
                              "Sig.Hypo.overlap") 
cat("\n")

"Total numbers of Promoters"
length(prom.gr)
cat("\n")

"Total numbers of Tiles"
length(myDiff.tiles.gr)
cat("\n")

"numbers of Tiles overlapping with a Promoter -1500 +1000"
sapply(difmeth.inter.PROM, nrow)
cat("\n")

"% of Tiles overlapping with a Promoter -1500 +1000"
round(100 * sapply(difmeth.inter.PROM, nrow)/nrow(myDiff.tiles) , 2)
cat("\n")

"Head of Tiles overlapping with a Promoter -1500 +1000 (GR List of 4 Genomic Ranges)"
lapply(difmeth.inter.PROM, head)
cat("\n")

"Writting MethylKitDMT_all.overlap_PROM "
for(i in 1:length(difmeth.inter.PROM) ) {
  write.table(difmeth.inter.PROM[[i]],
              paste0(REDIR, "MethylKitDMT_",
                     names(difmeth.inter.PROM[i]),
                     "_Promoters_myDiff_delta", 
                     DIF, "q", QV*100,"pc_", CONDI,".tsv"),
              sep='\t',
              dec=",",
              col.names=T,
              row.names=F,
              quote=F) 
}

cat("\n")
"Check the 4 outputs as ..."
"MethylKitDMT_all.overlap_Promoters_myDiff_xxxxxxxxxx.tsv"
"MethylKitDMT_Sig.overlap_Promoters_myDiff_xxxxxxxxxx.tsv"
"MethylKitDMT_HyperSig.overlap_Promoters_myDiff_xxxxxxxxxx.tsv"
"MethylKitDMT_HypoSig.overlap_Promoters_myDiff_xxxxxxxxxx.tsv"
cat("\n")
cat("\n")

"#################################################################"
"Extraction & Annotation of myDiff Tiles that overlap with Gene Bodies"
"#################################################################"
difmeth = myDiff.tiles.grl
annot = mm10Genes.gr
type = "GeneBody"

difmeth.inter.GB <- lapply(difmeth, ExportDiffAnnot, annot, type)

names(difmeth.inter.GB) <- c("all.overlap",
                              "Sig.overlap",
                              "Sig.Hyper.overlap",
                              "Sig.Hypo.overlap") 
cat("\n")

"Total numbers of Genes"
length(mm10Genes.gr)
cat("\n")

"Total numbers of Tiles"
length(myDiff.tiles.gr)
cat("\n")

"numbers of Tiles overlapping with a GeneBody"
sapply(difmeth.inter.GB, nrow)
cat("\n")

"% of Tiles overlapping with a GeneBody"
round(100 * sapply(difmeth.inter.GB, nrow)/nrow(myDiff.tiles) , 2)
cat("\n")

"Head of Tiles overlapping with a GeneBody (GR List of 4 Genomic Ranges)"
lapply(difmeth.inter.GB, head)
cat("\n")

"Writting MethylKitDMT_all.overlap_GeneBody "
for(i in 1:length(difmeth.inter.GB) ) {
  write.table(difmeth.inter.GB[[i]],
              paste0(REDIR, "MethylKitDMT_",
                     names(difmeth.inter.GB[i]),
                     "_GeneBody_myDiff_delta", 
                     DIF, "q", QV*100,"pc_", CONDI,".tsv"),
              sep='\t',
              dec=",",
              col.names=T,
              row.names=F,
              quote=F) 
}

cat("\n")
"Check the 4 outputs as ..."
"MethylKitDMT_all.overlap_GeneBody_myDiff_xxxxxxxxxx.tsv"
"MethylKitDMT_Sig.overlap_GeneBody_myDiff_xxxxxxxxxx.tsv"
"MethylKitDMT_HyperSig.overlap_GeneBody_myDiff_xxxxxxxxxx.tsv"
"MethylKitDMT_HypoSig.overlap_GeneBody_myDiff_xxxxxxxxxx.tsv"
cat("\n")
cat("\n")



cat("\n")
"#################################################################"
cat("\n")

" report SessionInfo"
sessionInfo()

cat("\n")
"#################################################################"
cat("\n")

# save the R Session in RData folder
save.image(file =paste0("../RData/", "MethylKitDMT_", CONDI,"_", TS,"_", RTS,"_", COV,"_delta", DIF,"_q", QV*100,"pc.RData") ) 
"An RData image is save in ../RData folder"


cat("\n")
"#################################################################"
" The END"
"#################################################################"
cat("\n")
"Here it's done. If you can read this it means the methylkit DMTiles script ran properly"
" for each differential analysis, 1 .jpg, 3 .pdf, 5 .bed, 12 .tsv (annotations) and 2 or 3 .xls(or txt) files are generated"
cat("\n")
cat("\n")
" - 1 output file with some statistics and parameters reminders"
cat("\n")
cat("\n")
" Files for Visualization"
cat("\n")
" 1 JPG file"
" - MethylKitDMT_VolcanoPlot_MethDiff_CondvsCtrl.jpg"
cat("\n")
" 3 PDF files"
" - MethylKitDMT_Pvalue_Distribution_MethDiff_CondvsCtrl.pdf"
" - MethylKitDMT_Qvalue_Distribution_MethDiff_CondvsCtrl.pdf"
" - MethylKitDMT_ChromPlot_MethDiff_CondvsCtrl.pdf"
cat("\n")
" 5 bed files"
" - MethylKitDMT_myDiff.tiles_CondvsCtrlTileSizeCov.bed "
" - MethylKitDMT_CondvsCtrl_mydeltaXXqXpc.bed"   
" - MethylKitDMT_CondvsCtrl_delta25q5pc.bed "
" - MethylKitDMT_CondvsCtrl_delta25q1pc.bed "
" - MethylKitDMT_CondvsCtrl_delta25q01pc.bed "
cat("\n")
cat("\n")
" Files for analysis, report and exploration"
cat("\n")
" 2 or 3 xls files"
" - MethylKitDMT_all_CondvsCtrl.txt.xls (with all tiles, txt file)"
" - MethylKitDMT_signif_CondvsCtrl.xls (with significant tiles, 3 sheets) or MethylKitDMT_signif_Diff25q01pc_CondvsCtrl.txt.xls only, if more than 65530 lines (txt file)"
" - MethylKitDMT_myDiff_deltaXXqXXpc_CondvsCtrl.txt.xls if more than 65530 lines (txt file)"
cat("\n")
" 12 annotations tsv files"
"  4 outputs as with CGI"
" - MethylKitDMT_all.overlap_CGI_myDiff_xxxxxxxxxx.tsv"
" - MethylKitDMT_Sig.overlap_CGI_myDiff_xxxxxxxxxx.tsv"
" - MethylKitDMT_HyperSig.overlap_CGI_myDiff_xxxxxxxxxx.tsv"
" - MethylKitDMT_HypoSig.overlap_CGI_myDiff_xxxxxxxxxx.tsv"
"  4 outputs with Promoters"
" - MethylKitDMT_all.overlap_Promoters_myDiff_xxxxxxxxxx.tsv"
" - MethylKitDMT_Sig.overlap_Promoters_myDiff_xxxxxxxxxx.tsv"
" - MethylKitDMT_HyperSig.overlap_Promoters_myDiff_xxxxxxxxxx.tsv"
" - MethylKitDMT_HypoSig.overlap_Promoters_myDiff_xxxxxxxxxx.tsv"
"  4 outputs with GeneBody"
" - MethylKitDMT_all.overlap_GeneBody_myDiff_xxxxxxxxxx.tsv"
" - MethylKitDMT_Sig.overlap_GeneBody_myDiff_xxxxxxxxxx.tsv"
" - MethylKitDMT_HyperSig.overlap_GeneBody_myDiff_xxxxxxxxxx.tsv"
" - MethylKitDMT_HypoSig.overlap_GeneBody_myDiff_xxxxxxxxxx.tsv"


cat("\n")
"Enjoy all these amazing results!!"

cat("\n")
"######################" 
cat("\n")

"Start"
start.time
"End"
Sys.time()

# Close output file
sink()


