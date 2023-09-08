---
title: "Code for the analysis in R"
author: "zhangyl"
date: '2023-09-07'
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

** this is the code for the analysis of single cell ATAC sequencing in the paper
create seurat object, quality control and visualition of data in 1k-cell stage

## R Markdown
```{r}
library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(ggplot2)
library(patchwork)
library(hdf5r)
library(dplyr)
library(readr)
library(pheatmap)
library(ggrepel)
library(LSD)
library(MASS)

mtx_dir_path <- "~/zbf_scATAC_mtx_2/1k"
mtx_path <- paste(mtx_dir_path, "count_matrix_over_aggregate.mtx" , sep = '/')
feature_path <- paste(mtx_dir_path, "count_matrix_over_aggregate.rows" , sep = "/")
barcode_path <- paste(mtx_dir_path, "count_matrix_over_aggregate.cols" , sep = '/')
metadata_path <- paste(mtx_dir_path , "sample_info.csv" , sep = '/')

features <- readr::read_tsv(feature_path, col_names = F) %>% tidyr::unite(feature)
barcodes <- readr::read_tsv(barcode_path, col_names = F) %>% tidyr::unite(barcode)

mtx <- Matrix::readMM(mtx_path) %>%
  magrittr::set_rownames(features$feature) %>%
  magrittr::set_colnames(barcodes$barcode)

chrom_assay <- CreateChromatinAssay(
  counts = mtx,
  sep = c("_" , "_"),
  genome = 'danRer11',
  fragments = paste(mtx_dir_path , 'aggregate_fragments.tsv.gz' , sep = '/'),
)

atac_1k <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  project = 'zbfish_scATAC_1k',
  meta.data = metadata
)

# Annotation
ensdb.zbf <- ensDbFromGtf(gtf = "~/danRer11_annotation/Danio_rerio.GRCz11.108.gtf.gz")
ensdb.zbf <- EnsDb(ensdb.zbf)
annotation.zbf <- GetGRangesFromEnsDb(ensdb.zbf)
seqlevelsStyle(annotation.zbf) <- 'UCSC'
genome(annotation.zbf) <- "danRer11"

Annotation(atac_1k) <- annotation.zbf

## add plate
plate <- data.frame("p1")
i = 1
for(i in 1:ncol(atac_1k))
{
  s <- strsplit(colnames(atac_1k) , split = "_")[[i]]
  plate[i] <- paste(s[2] , "_" ,
                    s[3] ,sep = "")
}
atac_1k$plate <- factor(plate , levels = c("1k_p5","1k_p6","1k_p7","1k_p8","1k_p13","1k_p14"))
 
## add batch
atac_1k$batch <- as.character(atac_1k$plate)
atac_1k$batch[atac_1k$batch %in% c('1k_p5' , '1k_p6' , '1k_p7' , '1k_p8')] = "1k_batch1"
atac_1k$batch[atac_1k$batch %in% c('1k_p13' , '1k_p14')] = "1k_batch2"
##### Quality Control (QC)
# compute nucleosome signal per cell
atac_1k <- NucleosomeSignal(object = atac_1k)
# compute TSS enrichment score per cell
atac_1k <- TSSEnrichment(object = atac_1k, fast = FALSE)

VlnPlot(atac_1k , features = c("nFeature_peaks" , "nCount_peaks" , 
                                 "mt_content" , "mapping_rate" , 
                                 "sequencing_depth" , "uniq_nuc_frags" ,
                                 "nucleosome_signal" , "TSS.enrichment") ,ncol = 2, group.by = "batch")
control.list <- grep("_384" , rownames(atac_1k@meta.data) , value = T)
## density heatmap
# LSD
heatscatter(log10(atac_1k$uniq_nuc_frags) , atac_1k$TSS.enrichment)
# ggplot
x = log10(atac_1k$uniq_nuc_frags)
y = atac_1k$TSS.enrichment
dens <- kde2d(x , y)
gr <- data.frame(with(dens , expand.grid(x , y)) , as.vector(dens$z))
names(gr) <- c("xgr" , "ygr" , "zgr")
mod <- loess(data = gr , zgr~xgr*ygr)
pointdens <- predict(mod , newdata = data.frame(xgr = x , ygr = y))
bar <- names(pointdens)
d <- as.data.frame(cbind(x , y , pointdens))
d$bar <- bar
colnames(d) <- c("log10_unique_fragments" , "TSS_Enrichment" , "Density" , "Barcode")

ggplot(data = d , aes(x = log10_unique_fragments , y = TSS_Enrichment , color = Density))+
  geom_point()+
  geom_hline(yintercept = 1 , linetype = "dashed")+
  geom_vline(xintercept = 3 , linetype = "dashed")+
  theme_bw()+
  theme(legend.background = element_blank() , panel.grid.major = element_blank() , panel.grid.minor = element_blank())+
  geom_label_repel(data = d[control.list,] , aes(label = Barcode),fontface = "bold" , color = "black",size = 3 , segment.color = "black",
                   point.padding = 1 , box.padding = 1)+
  scale_colour_gradientn(colours = c("grey" , "grey" , rainbow(7)[5:1]))

#### set a low nucleosome_signal value for 1k stage because it has a great influence on batch effect. 
atac_1k$nucleosome_group <- ifelse(atac_1k$nucleosome_signal > 3 , 'NS > 3' , "NS < 3")
FragmentHistogram(atac_1k , group.by = "nucleosome_group")

atac_1k <- RunTFIDF(atac_1k)
atac_1k <- FindTopFeatures(atac_1k, min.cutoff = 'q0')
atac_1k <- RunSVD(atac_1k)

atac_1k <- RunUMAP(object = atac_1k, reduction = 'lsi', dims = 2:30)
atac_1k <- FindNeighbors(object = atac_1k, reduction = 'lsi', dims = 2:30)
atac_1k <- FindClusters(object = atac_1k, verbose = FALSE, algorithm = 3)
DimPlot(object = atac_1k, label = TRUE) + NoLegend()
DimPlot(object = atac_1k, label = TRUE , group.by = "batch")
DimPlot(object = atac_1k, label = TRUE , group.by = "plate")
FeaturePlot(atac_1k , features = c("nucleosome_signal" , "TSS.enrichment"))

######## remove cells that are outliers for these QC metrics.
atac_1k <- subset(atac_1k , 
                    subset = uniq_nuc_frags > 1000 & 
                      TSS.enrichment > 1 & 
                      nucleosome_signal < 3 &
                      mapping_rate > 70)

atac_1k <- RunTFIDF(atac_1k)
atac_1k <- FindTopFeatures(atac_1k, min.cutoff = 'q0')
atac_1k <- RunSVD(atac_1k)

atac_1k <- RunUMAP(object = atac_1k, reduction = 'lsi', dims = 2:50)
DimPlot(object = atac_1k, label = TRUE) + NoLegend()
DimPlot(object = atac_1k, label = TRUE , group.by = "batch")
DimPlot(object = atac_1k, label = TRUE , group.by = "plate")

saveRDS(object = atac_1k , file = "~/zbf_scATAC_n/atac_1k_qc.RDS")
```

This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:

```{r cars}
summary(cars)
```

## Including Plots

You can also embed plots, for example:

```{r pressure, echo=FALSE}
plot(pressure)
```

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
