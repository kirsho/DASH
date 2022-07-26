#!/usr/bin/env Rscript 

library(magrittr)
library(dplyr)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(annotatr)


setwd("C:/Users/olivi/Nextcloud/Bioinfo_Analyses/WGBSDROM_script/script/figures_V2/")
mydir <- getwd()

load("Annotatr_mm10.RData") 


#################################################################
# add.annot
#################################################################

add.annot <-  function(filepath, genome, trackname){
  
  mybed <-  filepath
  mygenome <-  genome
  myname <-  trackname
  mybeddf <- read.table(mybed, header = F)
  write.table( mybeddf[,1:3], "mybeddf", col.names = F, row.names = F, quote = F, sep = "\t")
  mytrack <- "mybeddf"
  custom <- paste0(genome,"_custom_", myname)
  
  read_annotations(con = mytrack, 
                   genome = mygenome ,
                   name = myname,
                   format = 'bed')
  
  rm( mybed,
      mygenome,
      myname,
      mybeddf,
      mytrack)
  
  file.remove("mybeddf")
  
}

#################################################################
# Custom tracks
#################################################################
print("Reading and adding Custom tracks")


#Cluster
path2file <- "C:/Users/olivi/Nextcloud/Bioinfo_Analyses/WGBSDROM_script/script/figures_V2/preps/Custom_tracks__V2/custom_tracks/"

setwd(path2file)



for(file in dir()){
  print(file)
  
  add.annot(
      filepath = file,
      genome = "mm10",
      trackname = gsub(".BED", "", file )
    
  )

}


#################################################################
## Edit Annots  & annotations GR


( Annots <- c(Annots, annotatr_cache$list_env() ) )

( annotations <- build_annotations( genome = ORG, annotations = Annots ) )




#################################################################
# Edit your Custom shortcuts
#################################################################
cat("default")
cat("\n")

(Default <- c(  paste0( ORG , "_genes_intergenic"),
                paste0( ORG , "_genes_promoters"),
                paste0( ORG , "_genes_exons"),
                paste0( ORG , "_genes_introns"),
                paste0( ORG , "_cpg_islands"),
                paste0( ORG , "_cpg_shores"),
                paste0( ORG , "_cpg_shelves"),
                paste0( ORG , "_enhancers_fantom")
                )
)
                
                
              
cat("\n")


cat("Regulatory")
cat("\n")

(Regulatory <- c(paste0( ORG , "_genes_intergenic"),
                 paste0( ORG , "_genes_promoters"),
                 paste0( ORG , "_enhancers_fantom"),
                 paste0( ORG , "_custom_Promotersupregulated"),
                 paste0( ORG , "_custom_Promotersdownregulated"),
                 paste0( ORG , "_custom_Promotersbivalent"),
                 paste0( ORG , "_custom_PromotersK4K27"),      
                 paste0( ORG , "_custom_PromotersBEND3Science"),   
                 paste0( ORG , "_custom_EnhancersBEND3Science"), 
                 paste0( ORG , "_custom_PromotersBEND3"),   
                 paste0( ORG , "_custom_EnhancersBEND3"), 
                 paste0( ORG , "_custom_H3K4me3Promoters"),   
                 paste0( ORG , "_custom_H3K27me3Promoters")
                )
)
                 

cat("\n")


cat("BEND3")
cat("\n")

(BEND3 <- c(     paste0( ORG , "_genes_intergenic"),
                 paste0( ORG , "_custom_BEND3Science"),   
                 paste0( ORG , "_custom_BEND3Science1kbflanking"), 
                 paste0( ORG , "_custom_BEND3"),   
                 paste0( ORG , "_custom_BEND31kbflanking"),
                 paste0( ORG , "_custom_H3K4me3BEND3Science"),   
                 paste0( ORG , "_custom_H3K4me3BEND3"), 
                 paste0( ORG , "_custom_H3K27me3BEND3Science"),   
                 paste0( ORG , "_custom_H3K27me3BEND3")
            )
)  
                 
cat("\n")


cat("Histones")
cat("\n")

(Histones <- c(  paste0( ORG , "_genes_intergenic"),
                 paste0( ORG , "_custom_H3K4me3"),   
                 paste0( ORG , "_custom_H3K4me3BEND3Science"),   
                 paste0( ORG , "_custom_H3K4me3BEND3"), 
                 paste0( ORG , "_custom_H3K4me3K27me3"), 
                 paste0( ORG , "_custom_H3K4me3Promoters"),  
                 paste0( ORG , "_custom_H3K27me3"),   
                 paste0( ORG , "_custom_H3K27me3BEND3Science"),   
                 paste0( ORG , "_custom_H3K27me3BEND3"),  
                 paste0( ORG , "_custom_H3K27me3K4me3"), 
                 paste0( ORG , "_custom_H3K27me3Promoters")
              )
)  

cat("\n")



## All Customs features ( do not edit)
cat("Customs")
cat("\n")

(Customs <- grep("_custom_", Annots , value = T))

#################################################################
## Edit Update shortcuts            
shortcuts <- c(shortcuts,                    
               list(Customs=Customs,
                    Default=Default,       ## edit
                    Regulatory=Regulatory, ## edit
                    BEND3=BEND3,           ## edit
                    Histones=Histones      ## edit
               ) )


#################################################################
## Save image 
setwd(mydir)

save.image("Annotatr_mm10_Customs.RData")
