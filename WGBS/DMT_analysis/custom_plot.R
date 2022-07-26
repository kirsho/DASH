#!/usr/bin/env Rscript

#  Custom_plots for Lounis Yakhou Drom projects
# v1.1 2022-07-22 for 2i
## Initialise output file
sink(file = paste0(REDIR, CONDI, "_custom_plots_output.txt") )

# run in ifb from previous Figure RData
#rm(
#bpst,
#piedm,
#piedms,
#bpct,
#bpcd,
#bpsdiff,
#sas,
#saa,
#histmeth,
#histmeths,
#vp2,
#vps2,
#vc,
#dm_annotated,
#pcf,
#pcs
# )

cat("\n")
"#################################################################"
"Some Stats" %>% print
"#################################################################"
"Ana" %>% print
ANA %>% print
"Deciles Ctrl"  %>% print
AllDiffTiles$Methyl.Ctrl %>% quantile(., probs = seq(0, 1, 1/10) ) %>% print
cat("\n") 
"Quartiles Ctrl"  %>% print
AllDiffTiles$Methyl.Ctrl %>% quantile(., probs = seq(0, 1, 1/4) ) %>% print
cat("\n") 
"Deciles Conditions"  %>% print
AllDiffTiles$Methyl.Cond %>% quantile(., probs = seq(0, 1, 1/10) ) %>% print
cat("\n") 
"Quartiles Conditions" %>% print
AllDiffTiles$Methyl.Cond %>% quantile(., probs = seq(0, 1, 1/4) ) %>% print


"Counts Methstatus in AllDiffTiles" %>% print
"Ctrl"
AllDiffTiles$Methstatus_Ct %>% table %>% print
cat("\n")
AllDiffTiles$Methstatus_Ct  %>% table / (AllDiffTiles %>% nrow)  %>% print
cat("\n")
"Conditions"
AllDiffTiles$Methstatus_Cd %>% table %>% print
cat("\n")
AllDiffTiles$Methstatus_Cd  %>% table / (AllDiffTiles %>% nrow) %>% print

cat("\n")
"#################################################################"
"Plots"
"#################################################################"
cat("\n")

"#################################################################"
"Barplot Methylation Status bpst object"
"requires methst object"

bpst <- ggplot(data=methst , 
                aes( x = ID, fill = factor( Methstatus, rev( levels(Methstatus) ) ) )
                ) +  
        geom_bar() +
        labs(fill = '') +
        scale_fill_grey() +
        ggtitle(paste0(ANAtitle, " Methylation status") ) +
        theme_bw()  +
        theme(axis.title.x=element_blank(),
              plot.title = element_text(hjust = 0.5))

if(EXTR == "yes"){
svg( paste0(REDIR, ANA, "_Methylation_Status_barplot.svg") ) 
    print(bpst)
dev.off()
                }
if(EXTR == "no")
print("object saved in rdata")

"#################################################################"
"Riverplot Methylation Status rvpt object"

dd <- dfpc[, c("ID","Methstatus_Ct", "Methstatus_Cd")] %>% table %>% transform
dd$Methstatus_Ct <- factor(dd$Methstatus_Ct,  rev( levels(dd$Methstatus_Ct)))
dd$Methstatus_Cd <- factor(dd$Methstatus_Cd,  rev( levels(dd$Methstatus_Cd)))

rvpt <- ggplot(data= dd,
            aes( y = Freq,
                 axis1 = Methstatus_Ct,
                 axis2 = Methstatus_Cd,
                 label = ID ) ) +

      geom_stratum(size=1.2) +
      #geom_lode(aes(fill =  Methstatus_Ct), show.legend = T) +
      geom_alluvium(aes(fill = Methstatus_Ct),
                    knot.pos = 1/6) +   
      scale_fill_grey(start=0,end =0.85) +
      geom_text(color= "darkorange",fontface=2, stat = "stratum", aes(label = after_stat(stratum))) +
      scale_x_discrete(limits = sample.ids ) +
      labs(fill = 'Methylation status',
           y = "Counts" ) +
      ggtitle(paste0(ANAtitle, " Riverplot") ) +
      theme_light() +
      theme(axis.title.x=element_blank(),
            plot.title = element_text(hjust = 0))  


if(EXTR == "yes"){
svg( paste0(REDIR, ANA, "_Methylation_Status_Riverplot.svg") ) 
    print(rvpt)
dev.off()
                }
if(EXTR == "no")
print("object saved in rdata")



"#################################################################"
"Riverplot Methylation Status on annots"
"requires methst object"

IGV = c("mm10_genes_intergenic", "mm10_genes_promoters", "mm10_enhancers_fantom", 
		"mm10_cpg_islands", "mm10_cpg_shores", "mm10_cpg_shelves", 
		"mm10_custom_H3K4me3", "mm10_custom_H3K27me3", "mm10_custom_BEND3")

for(i in IGV ){
  print(i)

  dfpc_temp <- dfpc_annotated %>% data.frame
  dfpc_temp <-  dfpc_temp[ which( dfpc_temp$annot.type %in% i ) , ]
  dfpc_temp$annot.type <- factor(dfpc_temp$annot.type, levels = IGV[i] )

  dd <- dfpc_temp[, c("ID","Methstatus_Ct", "Methstatus_Cd")] %>% table %>% transform
  dd$Methstatus_Ct <- factor(dd$Methstatus_Ct,  rev( levels(dd$Methstatus_Ct)))
  dd$Methstatus_Cd <- factor(dd$Methstatus_Cd,  rev( levels(dd$Methstatus_Cd)))

  if(EXTR == "yes"){

                    rvpt2 <- ggplot(data= dd,
                              aes( y = Freq,
                                   axis1 = Methstatus_Ct,
                                   axis2 = Methstatus_Cd,
                                   label = ID ) ) +

                              geom_stratum(size=1.2) +
                              #geom_lode(aes(fill =  Methstatus_Ct), show.legend = T) +
                              geom_alluvium(aes(fill = Methstatus_Ct),
                                            knot.pos = 1/6) +   
                              scale_fill_grey(start=0,end =0.85) +
                              geom_text(color= "darkorange",fontface=2, stat = "stratum", aes(label = after_stat(stratum))) +
                              scale_x_discrete(limits = sample.ids ) +
                              labs(fill = 'Methylation status',
                                   y = "Counts" ) +
                              ggtitle(paste0(ANAtitle, " Riverplot"),
                                      subtitle = paste0("On annotation track ", i) ) +
                              theme_light() +
                              theme(axis.title.x=element_blank(),
                                    plot.title = element_text(hjust = 0))  
                    

                    svg( paste0(REDIR, ANA, "_MethStatus_by_Methlevel_Riverplot_in_", i , ".svg") )
                    print( rvpt2 )
                    dev.off()

                    rm( rvpt2, dfpc_temp, dd )

                }
  
  if(EXTR == "no"){
    print("object saved in rdata")
                  }
  
}





"#################################################################"
"Piechart Differentially methylated Tiles distribution "
"piedm object on AllDiffTiles$DM_status %>% table %>% data.frame "
piedm <- ggplot(AllDiffTiles$DM_status %>% table %>% data.frame,
                aes(x="", y=Freq, fill=.)) +
          geom_col() +      
          scale_fill_manual("", values = c( 'Hyper' = "firebrick3",
                                            'Hypo'  = "chartreuse3",
                                            'None'  = "azure3") 
                            ) +
          coord_polar("y", start=0) +
          #geom_text(aes(label = Freq),
          #          position = position_stack(vjust = 0.5), col = "white", fontface = "bold") +
          theme_void() +
          ggtitle(paste0(ANAtitle),
                  subtitle = paste0(" Differentially methylated tiles distribution (all)")
                  ) +
          theme(plot.title = element_text(hjust = 0.5))

"piedms object on  filter(AllDiffTiles, DM_status != \"None\")  -> cnt "

filter(AllDiffTiles, DM_status != "None") %>% dplyr::select(., DM_status) %>% table %>% data.frame -> cnt

piedms <- ggplot( cnt[-3,], aes(x="", y=Freq, fill=.) ) +
          geom_col() +      
          scale_fill_manual("", values = c( 'Hyper' = "firebrick3",
                                            'Hypo'  = "chartreuse3") 
                            ) +
          coord_polar("y", start=0) +
          #geom_text(aes(label = Freq),
          #          position = position_stack(vjust = 0.5), col = "white", fontface = "bold") +
          theme_void() +
          ggtitle(paste0(ANAtitle),
                  subtitle = paste0(" Differentially methylated tiles distribution (Significant)")
                  ) +
          theme(plot.title = element_text(hjust = 0.5))


if(EXTR == "yes"){
svg( paste0(REDIR, ANA, "_DiffMeth_Tiles_distribution_piechart.svg") ) 
    grid.arrange(piedm, piedms, ncol = 1)
dev.off()
                  }
if(EXTR == "no")
print("object saved in rdata")

#
#"#################################################################"
#"barplot fraction low mid high in Diff meth"
#"bpct"
#bpct <- ggplot(data=AllDiffTiles , 
#              aes(x=Methstatus_Ct, fill=DM_status)) +
#        geom_bar() +
#        theme_bw() +
#        ggtitle(paste0(ANAtitle, " All-Differential Methylation % ~ Tiles methylation status in ", sample.ids[1]) ) +
#        theme(plot.title = element_text(hjust = 0.5)) 
#
#"bpct"
#bpcd <- ggplot(data=AllDiffTiles , 
#              aes(x=Methstatus_Cd, fill=DM_status)) +
#        geom_bar() +
#        theme_bw() +
#        ggtitle(paste0(ANAtitle, " All Differential Methylation % ~ Tiles methylation status in ", sample.ids[2]) ) +
#        theme(plot.title = element_text(hjust = 0.5)) 
#
#"bpsdiff"
#bpsdiff <- ggplot(data=filter(AllDiffTiles, Diff_expr =="Significant") , 
#              aes(x=Methstatus_Ct, fill=DM_status)) +
#        geom_bar() +
#        theme_bw() +
#        ggtitle(paste0(ANAtitle, " Signif Differential Methylation % ~ Tiles methylation status in ", sample.ids[1]) ) +
#        theme(plot.title = element_text(hjust = 0.5))      
#
#if(EXTR == "yes"){
#svg( paste0(REDIR, ANA, "_Fraction_MethStatus_in_DiffMeth_barplot.svgg"))  
#    grid.arrange(bpct, bpcd, bpsdiff,  ncol = 1)
#dev.off()
#                }
#
#if(EXTR == "no"){
#print("object saved in rdata")
#                }


"#################################################################"
" New barplot fraction low mid high in Diff meth on All the Annots"

for(i in IGV ){
  print(i)

  dfpc_temp <- dfpc_annotated %>% data.frame
  dfpc_temp <-  dfpc_temp[ which( dfpc_temp$annot.type %in% i ) , ]
  dfpc_temp$annot.type <- factor(dfpc_temp$annot.type, levels = IGV[i] )

  if(EXTR == "yes"){

                    bp2ct <-  ggplot( dfpc_temp, aes(x=ID, fill = DM_status)) +
                              geom_bar( ) +
                              facet_grid(~ Methstatus) +
                              theme_bw() +
                              scale_fill_manual("", values = c( 'Hyper' = "firebrick3",
                                                                'Hypo'  = "chartreuse3",
                                                                'None'  = "azure3") 
                                                ) +
                              ggtitle(paste0(label = ANAtitle, " Differential methylation status"),
                                             subtitle = paste0( "On annotation track ", i, 
                                                                " and by tiles methylation level in ",
                                                                sample.ids[1])
                                      ) +
                              theme(
                                axis.title.x=element_blank(),
                                plot.title = element_text(hjust = 0.5),
                                plot.subtitle = element_text(hjust = 0.5),
                                strip.text.x = element_text(size = 12, color = "black", face = "bold.italic")
                                    )

                    svg( paste0(REDIR, ANA, "_MethStatus_by_Methlevel_barplot_in_", i , ".svg") )
                    print( bp2ct )
                    dev.off()

                    rm( bp2ct, dfpc_temp )

                }
  
  if(EXTR == "no"){
    print("object saved in rdata")
                  }
  
}





"#################################################################"
"# Old barplot fraction low mid high in Diff meth on myshortcut BEND3"
#
# for(i in BEND3 ){
# print(i)
#
#
# Adtdf <- AllDiffTiles_annotated %>% data.frame
# Adtdf <- Adtdf[ which( Adtdf$annot.type %in% i ) , ]
# 
# x_order  <-  c(
#   'Hyper',
#   'Hypo' ,
#   'None')
# 
# if(EXTR == "yes"){
# 
# "bpctsct"
# bpctsct <- ggplot(data=Adtdf , 
#                   aes(x=Methstatus_Ct, fill=DM_status)) +
#   geom_bar() +
#   theme_bw() +
#   ggtitle(paste0(ANAtitle, " All-Differential Methylation % ~ Tiles methylation status in ", sample.ids[1], " on ", i ) ) +
#   theme(plot.title = element_text(hjust = 0.5)) 
# 
# "bpcdsct"
# bpcdsct <- ggplot(data=Adtdf , 
#                   aes(x=Methstatus_Cd, fill=DM_status)) +
#   geom_bar() +
#   theme_bw() +
#   ggtitle(paste0(ANAtitle, " All Differential Methylation % ~ Tiles methylation status in ", sample.ids[2], " on ", i ) ) +
#   theme(plot.title = element_text(hjust = 0.5)) 
# 
# "bpsdiffsct"
# bpsdiffsct <- ggplot(data=filter(Adtdf, Diff_expr =="Significant") , 
#                      aes(x=Methstatus_Ct, fill=DM_status)) +
#   geom_bar() +
#   theme_bw() +
#   ggtitle(paste0(ANAtitle, " Signif Differential Methylation % ~ Tiles methylation status in ", sample.ids[1], " on ", i ) ) +
#   theme(plot.title = element_text(hjust = 0.5))      
# 
#
#   png( paste0(REDIR, ANA, "_Fraction_MethStatus_in_DiffMeth_barplot in ", i,  ".png"),
#        width = 1200, height = 1200, units = "px", pointsize = 12 )  
#   
#   print( grid.arrange(bpctsct, bpcdsct, bpsdiffsct,  ncol = 1) )
#   
#   dev.off()
#   
#   
#   rm( bpctsct, bpcdsct, bpsdiffsct, Adtf )
# }
# 
# if(EXTR == "no"){
#   print("object saved in rdata")
# }
# 
# }
#
#

"##################################################################"
"Count annotations" %>% print

summarize_annotations(
    annotated_regions = SigDiffTiles_annotated,
    annotated_random = random_SigDiffTiles_annotated) %>% data.frame -> sas

sas$annot.type <- factor(sas$annot.type, levels=Annots)

summarize_annotations(
    annotated_regions = AllDiffTiles_annotated,
    annotated_random = random_AllDiffTiles_annotated) %>% data.frame -> saa

saa$annot.type <- factor(saa$annot.type, levels=Annots)

cat("\n") 
"SigDiffTiles_annotated summarize_annotations" %>% print
sas %>% print

cat("\n") 
"AllDiffTiles_annotated summarize_annotations" %>% print
saa %>% print


"#################################################################"
" Histogram of methylation all"
"histmeth"
histmeth <- ggplot(dfpc %>% data.frame, 
                    aes(x=Meth_perc*100, fill=ID, color=ID)
                    ) +
            scale_fill_manual(values=c("darkorange2", "cornflowerblue")) +      
            geom_histogram(position="identity",
                           alpha=0.6,
                           breaks = seq(0, 100, 2)) +                           
            ggtitle(paste0("Methylation distribution on tiles (All)") ) +
            xlab("% of Methylation") + ylab("Counts") +
            theme_bw() +  
            theme(legend.title = element_blank(),
                  plot.title = element_text(hjust = 0.5)) 

if(EXTR == "yes"){
svg( paste0(REDIR, ANA, "_Histogram_of_Methylation_All.svg") )  
   print(histmeth) 
dev.off()
                }

if(EXTR == "no"){
print("object saved in rdata")
                }



"#################################################################"
" Histogram of methylation Signif "
" histmeths "
histmeths <-  ggplot(dfpc_Sig %>% data.frame, 
                     aes(x=Meth_perc*100, fill=ID, color=ID)
                    ) +             
              scale_fill_manual(values=c("darkorange2", "cornflowerblue")) +     
              geom_histogram(position="identity",
                             alpha=0.6,
                             breaks = seq(0, 100, 2)) +                           
              ggtitle(paste0("Methylation distribution on tiles (Significant)") ) +
              xlab("% of Methylation") + ylab("Counts") +
              theme_bw() +  
              theme(legend.title = element_blank(),
                    plot.title = element_text(hjust = 0.5))


if(EXTR == "yes"){
svg( paste0(REDIR, ANA, "_HistogramofMethylation_Signif.svg") )
   print(histmeths) 
dev.off()
                }

if(EXTR == "no"){
print("object saved in rdata")
                }
# before as png
#png( paste0(REDIR, ANA, "_HistogramofMethylation_Signif.png"),
#        width = 1200, height = 1200, units = "px", pointsize = 12 )                  


"#################################################################"
" Violin plot of methylation All "
" vp2 "
vp2 <- ggplot(dfpc, aes(x=ID, y=Meth_perc*100, fill = ID)) +  
        geom_violin(trim=FALSE) +  
        geom_boxplot(width=0.05) + 
        labs(title = " Tiles (all) methylation % distribution" ,
              x = "Samples", y = "Methylation %" )  +
        theme_bw() +               
        scale_fill_manual(values=c("darkorange2", "cornflowerblue")) +  
        theme(axis.title.x=element_blank(),
              plot.title = element_text(hjust = 0.5), 
              legend.position = "none")

if(EXTR == "yes"){
svg( paste0(REDIR, ANA, "_Violinplot_Methylation_All.svg") )  
   print(vp2) 
dev.off()
                }

if(EXTR == "no"){
print("object saved in rdata")
                }


"Violin plot of methylation Signif"
" vp2s "
vps2 <- ggplot(dfpc_Sig, aes(x=ID, y=Meth_perc*100, fill = ID)) +  
        geom_violin(trim=FALSE) +  
        geom_boxplot(width=0.05) + 
        labs( title = "Tiles (Significant) methylation % distribution" ,
              x = "Samples", y = "Methylation %" )  +
        theme_bw() +               
        scale_fill_manual(values=c("darkorange2", "cornflowerblue")) +  
        theme(axis.title.x=element_blank(),
              plot.title = element_text(hjust = 0.5), 
              legend.position = "none")


if(EXTR == "yes"){
svg( paste0(REDIR, ANA, "_Violinplot_Methylation_Signif.svg") ) 
   print(vps2) 
dev.off()
                }

if(EXTR == "no"){
print("object saved in rdata")
                }


"#################################################################"
" Volcano plot"
" vc "
vc <- ggplot(data = AllDiffTiles, 
              aes(x = meth.diff, y = -log10(qvalue), col = DM_status)) +
        geom_point() + 
        theme_bw() + 
        labs(col = '') +
        xlab("Differential methylation") +
        xlim(-100, 100) +
        scale_color_manual(values=c('Hyper' = "firebrick3",
                                    'Hypo' = "chartreuse3",
                                    'None' = "azure3")) +
        geom_vline(xintercept = c(-DIF, DIF), linetype="dashed", col="black") +
        geom_hline(yintercept = -log10(QV), linetype="dashed", col="black") +
        ggtitle(paste0(ANAtitle, " Volcano plot on Tiles")) +
        theme(plot.title = element_text(hjust = 0.5))


if(EXTR == "yes"){
svg( paste0(REDIR, ANA, "_VolcanoPlot.svg") )  
   print(vc) 
dev.off()
                }

if(EXTR == "no"){
print("object saved in rdata")
                }



### plots all volcanos
for(i in IGV ){
  print(i)

  Adtdf <- AllDiffTiles_annotated %>% data.frame
  Adtdf <- Adtdf[ which( Adtdf$annot.type %in% i ) , ]
  Adtdf$annot.type <- factor(Adtdf$annot.type, levels = IGV[i] )



  if(EXTR == "yes"){

      vc_temp <- ggplot(data = Adtdf, 
                    aes(x = meth.diff, y = -log10(qvalue), col = DM_status)) +
              geom_point() + 
              theme_bw() + 
              labs(col = '') +
              xlab("Differential methylation") +
              xlim(-100, 100) +
              scale_color_manual(values=c('Hyper' = "firebrick3",
                                          'Hypo' = "chartreuse3",
                                          'None' = "azure3")) +
              geom_vline(xintercept = c(-DIF, DIF), linetype="dashed", col="black") +
              geom_hline(yintercept = -log10(QV), linetype="dashed", col="black") +
              ggtitle(paste0(ANAtitle, " Volcano plot on Tiles"),
                      subtitle = paste0("On annotation track ", i ) ) +
              theme(plot.title = element_text(hjust = 0.5))


if(EXTR == "yes"){
svg( paste0(REDIR, ANA, "_VolcanoPlot_on_track_", i,".svg") )  
   print(vc_temp) 
dev.off()
                }

if(EXTR == "no"){
print("object saved in rdata")
                }

                }
   }








"#################################################################"
" plot categorical"
" pcf & pcs "

dm_annotated <- AllDiffTiles_annotated  %>% data.frame
dm_annotated$annot.type <- factor(dm_annotated$annot.type , levels = Annots )


"relative counts"
pcr <- ggplot(data=dm_annotated, aes(x=annot.type, fill=DM_status)) +
        geom_bar(position = "fill") +
        scale_fill_manual("",
                    values = c( 'Hyper' = "firebrick3",
                                'Hypo'  = "chartreuse3",
                                'None'  = "azure3") 
                          ) +
        scale_x_discrete("", labels = gsub("*mm10_", "", Annots) ) +
        theme_bw() +
        ggtitle( paste0('DM Status on annotation tracks (relative counts)'),
        subtitle = paste0(ANAtitle) ) +
        theme(axis.text.x = element_text( size = 8, angle = 45, vjust = 1, hjust = 1),
              axis.title.x = element_blank(), legend.title = element_blank()
              ) 


if(EXTR == "yes"){
svg( paste0(REDIR, ANA, "_plot_categorical_relative.svg") ) 
   print(pcr) 
dev.off()
                }

if(EXTR == "no"){
print("object saved in rdata")
                }


#x_order  <-  c(
#    'Hyper',
#    'Hypo',
#    'None')
#
#fill_order <- Annots

#pcf <- plot_categorical(
#    annotated_regions = dm_annotated,
#    fill='DM_status', x='annot.type',
#    fill_order = x_order, x_order = fill_order, position='fill',
#    plot_title =  paste0(ANAtitle, ' DM Status by annotation (relative counts)') ,
#    legend_title =  'Annotations',
#    x_label = '',
#    y_label = 'proportion')
#

"absolute counts"
pca <- ggplot(data=dm_annotated, aes(x=annot.type, fill=DM_status)) +
        geom_bar() +
        scale_fill_manual("",
                    values = c( 'Hyper' = "firebrick3",
                                'Hypo'  = "chartreuse3",
                                'None'  = "azure3") 
                          ) +
        scale_x_discrete("", labels = gsub("*mm10_", "", Annots) ) +
        theme_bw() +
        ggtitle( paste0('DM Status on annotation tracks (absolute counts)'),
        subtitle = paste0(ANAtitle) ) +
        theme(axis.text.x = element_text( size = 8, angle = 45, vjust = 1, hjust = 1),
              axis.title.x = element_blank(), legend.title = element_blank()
              ) 


if(EXTR == "yes"){
svg( paste0(REDIR, ANA, "_plot_categorical_counts.svg") ) 
   print(pca) 
dev.off()
                }
                
if(EXTR == "no"){
print("object saved in rdata")
                }



"#################################################################"
" plot categorical on shortcuts for BEND3 project (plot 7)"
" pcrsct & pcasct"

myshortcuts <- list(#AllAnnots = Annots,
                   Default = Default,
                   Regulatory = Regulatory, 
                   BEND3 = BEND3,
                   Histones = Histones     
                    )

for(i in names(myshortcuts) ){
    print(i)
    print(myshortcuts[[i]])  

    dm_annotated <- AllDiffTiles_annotated %>% data.frame
    dm_annotated <- dm_annotated[ which( dm_annotated$annot.type %in% myshortcuts[[i]] ) , ]
    dm_annotated$annot.type <- factor(dm_annotated$annot.type, levels = myshortcuts[[i]] )

 

    pcrsct <- ggplot(data=dm_annotated, aes(x=annot.type, fill=DM_status)) +
              geom_bar(position = "fill") +
              scale_fill_manual("",
                          values = c( 'Hyper' = "firebrick3",
                                      'Hypo'  = "chartreuse3",
                                      'None'  = "azure3") 
                                ) +
              scale_x_discrete("", labels = gsub("*mm10_", "", myshortcuts[[i]]) ) +
              theme_bw() +
              ggtitle( paste0( 'DM Status on ', i, ' annotation tracks (relative counts)'),
              subtitle = paste0(ANAtitle)
                      ) +
              theme(axis.text.x = element_text( size = 10, angle = 45, vjust = 1, hjust = 1),
                    axis.title.x = element_blank(), legend.title = element_blank()
                    ) 


    if(EXTR == "yes"){
    svg( paste0(REDIR, ANA, "_plot_categorical_relative_on_shortcuts_", i ,".svg") ) 
      print(pcrsct) 
    dev.off()
                    }

    pcasct <- ggplot(data=dm_annotated, aes(x=annot.type, fill=DM_status)) +
            geom_bar() +
            scale_fill_manual("",
                        values = c( 'Hyper' = "firebrick3",
                                    'Hypo'  = "chartreuse3",
                                    'None'  = "azure3") 
                              ) +
            scale_x_discrete("", labels = gsub("*mm10_", "", myshortcuts[[i]]) ) +
            theme_bw() +
            ggtitle( paste0( 'DM Status on ', i, ' annotation tracks (absolute counts)'),
            subtitle = paste0(ANAtitle)
                    ) +
            theme(axis.text.x = element_text( size = 10, angle = 45, vjust = 1, hjust = 1),
                  axis.title.x = element_blank(), legend.title = element_blank()
                  ) 

    if(EXTR == "yes"){
    svg( paste0(REDIR, ANA, "_plot_categorical_absolute_on_shortcuts_", i ,".svg") ) 
       print(pcasct) 
    dev.off()
                    }

    if(EXTR == "no"){
    print("object pcsct is not saved in R saved in rdata")
                    }

    rm(pcrsct, pcasct, dm_annotated)

}




"#################################################################"
" Write bed files "

SigDiffTiles_annotated %>% data.frame()

for(i in IGV ){
  print(i)

      df_temp <- SigDiffTiles_annotated %>% data.frame()
      df_temp <- df_temp[ which( df_temp$annot.type %in% i ) ,
                         c(   "seqnames",
                              "start",
                              "end"  ,            
                              "meth.diff",
                              "annot.id",
                              "annot.tx_id",
                              "annot.gene_id" ,
                              "annot.symbol")
                          ]

      write.table(df_temp,
            paste0(REDIR, CONDI, "_in_",i ,"_.bed" ), 
            sep='\t',
            dec=".",
            col.names=F,
            row.names=F,
            quote=F)

      rm(df_temp)

}









sessionInfo()

sink()


# x_order  <-  c(
#    'Hyper',
#    'Hypo' ,
#    'None')
#
#fill_order <- myshortcuts[[i]]
#
#pcsct <- plot_categorical(
#    annotated_regions = dm_annotated,
#    fill='DM_status', x='annot.type',
#    fill_order = x_order, x_order = fill_order, position='fill',
#    plot_title =  paste0(ANA, ' DM Status by Annotation relative Counts for ', i ) ,
#    legend_title = 'Annotations',
#    x_label = 'DM status',
#    y_label = 'proportion')  