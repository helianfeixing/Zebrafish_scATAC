---
title: "Code for the analysis in R"
author: "zhangyl"
date: '2023-09-07'
output: github_document
---



** this is the code for the analysis of single cell ATAC sequencing in the paper
create seurat object, quality control and visualition of data in 1k-cell stage

## R Markdown

```r
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
```

```
## Rows: 331218 Columns: 3
## ── Column specification ──────────────────────────────────────────────────────────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr (1): X1
## dbl (2): X2, X3
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

```r
barcodes <- readr::read_tsv(barcode_path, col_names = F) %>% tidyr::unite(barcode)
```

```
## Rows: 2304 Columns: 1
## ── Column specification ──────────────────────────────────────────────────────────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr (1): X1
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

```r
mtx <- Matrix::readMM(mtx_path) %>%
  magrittr::set_rownames(features$feature) %>%
  magrittr::set_colnames(barcodes$barcode)

chrom_assay <- CreateChromatinAssay(
  counts = mtx,
  sep = c("_" , "_"),
  genome = 'danRer11',
  fragments = paste(mtx_dir_path , 'aggregate_fragments.tsv.gz' , sep = '/'),
)
```

```
## Computing hash
```

```
## Error in CreateFragmentObject(path = fragments, cells = cells, validate.fragments = validate.fragments, : Not all cells requested could be found in the fragment file.
```

```r
atac_1k <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  project = 'zbfish_scATAC_1k',
  meta.data = metadata
)
```

```
## Error in CreateSeuratObject(counts = chrom_assay, assay = "peaks", project = "zbfish_scATAC_1k", : 找不到对象'chrom_assay'
```

```r
# Annotation
ensdb.zbf <- ensDbFromGtf(gtf = "~/danRer11_annotation/Danio_rerio.GRCz11.108.gtf.gz")
```

```
## Importing GTF file ... OK
## Processing genes ...
```

```
## Warning in ensDbFromGRanges(GTF, outfile = outfile, path = path, organism = organism, : I'm missing column(s): 'entrezid'. The
## corresponding database column(s) will be empty!
```

```
##  Attribute availability:
##   o gene_id ... OK
##   o gene_name ... OK
##   o entrezid ... Nope
##   o gene_biotype ... OK
## OK
## Processing transcripts ... 
##  Attribute availability:
##   o transcript_id ... OK
##   o gene_id ... OK
##   o transcript_biotype ... OK
##   o transcript_name ... OK
## OK
## Processing exons ... OK
## Processing chromosomes ... Fetch seqlengths from ensembl ... OK
## Processing metadata ... OK
## Generating index ... OK
##   -------------
## Verifying validity of the information in the database:
## Checking transcripts ... OK
## Checking exons ... OK
```

```r
ensdb.zbf <- EnsDb(ensdb.zbf)
annotation.zbf <- GetGRangesFromEnsDb(ensdb.zbf)
```

```
## Warning in .merge_two_Seqinfo_objects(x, y): The 2 combined objects have no sequence levels in common. (Use
##   suppressWarnings() to suppress this warning.)
```

```
## Warning in .merge_two_Seqinfo_objects(x, y): The 2 combined objects have no sequence levels in common. (Use
##   suppressWarnings() to suppress this warning.)

## Warning in .merge_two_Seqinfo_objects(x, y): The 2 combined objects have no sequence levels in common. (Use
##   suppressWarnings() to suppress this warning.)

## Warning in .merge_two_Seqinfo_objects(x, y): The 2 combined objects have no sequence levels in common. (Use
##   suppressWarnings() to suppress this warning.)

## Warning in .merge_two_Seqinfo_objects(x, y): The 2 combined objects have no sequence levels in common. (Use
##   suppressWarnings() to suppress this warning.)

## Warning in .merge_two_Seqinfo_objects(x, y): The 2 combined objects have no sequence levels in common. (Use
##   suppressWarnings() to suppress this warning.)

## Warning in .merge_two_Seqinfo_objects(x, y): The 2 combined objects have no sequence levels in common. (Use
##   suppressWarnings() to suppress this warning.)

## Warning in .merge_two_Seqinfo_objects(x, y): The 2 combined objects have no sequence levels in common. (Use
##   suppressWarnings() to suppress this warning.)

## Warning in .merge_two_Seqinfo_objects(x, y): The 2 combined objects have no sequence levels in common. (Use
##   suppressWarnings() to suppress this warning.)

## Warning in .merge_two_Seqinfo_objects(x, y): The 2 combined objects have no sequence levels in common. (Use
##   suppressWarnings() to suppress this warning.)

## Warning in .merge_two_Seqinfo_objects(x, y): The 2 combined objects have no sequence levels in common. (Use
##   suppressWarnings() to suppress this warning.)

## Warning in .merge_two_Seqinfo_objects(x, y): The 2 combined objects have no sequence levels in common. (Use
##   suppressWarnings() to suppress this warning.)

## Warning in .merge_two_Seqinfo_objects(x, y): The 2 combined objects have no sequence levels in common. (Use
##   suppressWarnings() to suppress this warning.)

## Warning in .merge_two_Seqinfo_objects(x, y): The 2 combined objects have no sequence levels in common. (Use
##   suppressWarnings() to suppress this warning.)

## Warning in .merge_two_Seqinfo_objects(x, y): The 2 combined objects have no sequence levels in common. (Use
##   suppressWarnings() to suppress this warning.)

## Warning in .merge_two_Seqinfo_objects(x, y): The 2 combined objects have no sequence levels in common. (Use
##   suppressWarnings() to suppress this warning.)

## Warning in .merge_two_Seqinfo_objects(x, y): The 2 combined objects have no sequence levels in common. (Use
##   suppressWarnings() to suppress this warning.)

## Warning in .merge_two_Seqinfo_objects(x, y): The 2 combined objects have no sequence levels in common. (Use
##   suppressWarnings() to suppress this warning.)

## Warning in .merge_two_Seqinfo_objects(x, y): The 2 combined objects have no sequence levels in common. (Use
##   suppressWarnings() to suppress this warning.)

## Warning in .merge_two_Seqinfo_objects(x, y): The 2 combined objects have no sequence levels in common. (Use
##   suppressWarnings() to suppress this warning.)

## Warning in .merge_two_Seqinfo_objects(x, y): The 2 combined objects have no sequence levels in common. (Use
##   suppressWarnings() to suppress this warning.)

## Warning in .merge_two_Seqinfo_objects(x, y): The 2 combined objects have no sequence levels in common. (Use
##   suppressWarnings() to suppress this warning.)

## Warning in .merge_two_Seqinfo_objects(x, y): The 2 combined objects have no sequence levels in common. (Use
##   suppressWarnings() to suppress this warning.)

## Warning in .merge_two_Seqinfo_objects(x, y): The 2 combined objects have no sequence levels in common. (Use
##   suppressWarnings() to suppress this warning.)

## Warning in .merge_two_Seqinfo_objects(x, y): The 2 combined objects have no sequence levels in common. (Use
##   suppressWarnings() to suppress this warning.)
```

```r
seqlevelsStyle(annotation.zbf) <- 'UCSC'
genome(annotation.zbf) <- "danRer11"

Annotation(atac_1k) <- annotation.zbf
```

```
## Error in Annotation(atac_1k) <- annotation.zbf: 找不到对象'atac_1k'
```

```r
## add plate
plate <- data.frame("p1")
i = 1
for(i in 1:ncol(atac_1k))
{
  s <- strsplit(colnames(atac_1k) , split = "_")[[i]]
  plate[i] <- paste(s[2] , "_" ,
                    s[3] ,sep = "")
}
```

```
## Error in h(simpleError(msg, call)): 在为'ncol'函数选择方法时评估'x'参数出了错: 找不到对象'atac_1k'
```

```r
atac_1k$plate <- factor(plate , levels = c("1k_p5","1k_p6","1k_p7","1k_p8","1k_p13","1k_p14"))
```

```
## Error in atac_1k$plate <- factor(plate, levels = c("1k_p5", "1k_p6", "1k_p7", : 找不到对象'atac_1k'
```

```r
## add batch
atac_1k$batch <- as.character(atac_1k$plate)
```

```
## Error in eval(expr, envir, enclos): 找不到对象'atac_1k'
```

```r
atac_1k$batch[atac_1k$batch %in% c('1k_p5' , '1k_p6' , '1k_p7' , '1k_p8')] = "1k_batch1"
```

```
## Error in atac_1k$batch[atac_1k$batch %in% c("1k_p5", "1k_p6", "1k_p7", : 找不到对象'atac_1k'
```

```r
atac_1k$batch[atac_1k$batch %in% c('1k_p13' , '1k_p14')] = "1k_batch2"
```

```
## Error in atac_1k$batch[atac_1k$batch %in% c("1k_p13", "1k_p14")] = "1k_batch2": 找不到对象'atac_1k'
```

```r
##### Quality Control (QC)
# compute nucleosome signal per cell
atac_1k <- NucleosomeSignal(object = atac_1k)
```

```
## Error in DefaultAssay(object = object): 找不到对象'atac_1k'
```

```r
# compute TSS enrichment score per cell
atac_1k <- TSSEnrichment(object = atac_1k, fast = FALSE)
```

```
## Error in DefaultAssay(object = object): 找不到对象'atac_1k'
```

```r
VlnPlot(atac_1k , features = c("nFeature_peaks" , "nCount_peaks" , 
                                 "mt_content" , "mapping_rate" , 
                                 "sequencing_depth" , "uniq_nuc_frags" ,
                                 "nucleosome_signal" , "TSS.enrichment") ,ncol = 2, group.by = "batch")
```

```
## Error in DefaultAssay(object = object): 找不到对象'atac_1k'
```

```r
control.list <- grep("_384" , rownames(atac_1k@meta.data) , value = T)
```

```
## Error in h(simpleError(msg, call)): 在为'grep'函数选择方法时评估'x'参数出了错: 在为'rownames'函数选择方法时评估'x'参数出了错: 找不到对象'atac_1k'
```

```r
## density heatmap
# LSD
heatscatter(log10(atac_1k$uniq_nuc_frags) , atac_1k$TSS.enrichment)
```

```
## Error in heatscatter(log10(atac_1k$uniq_nuc_frags), atac_1k$TSS.enrichment): 找不到对象'atac_1k'
```

```r
# ggplot
x = log10(atac_1k$uniq_nuc_frags)
```

```
## Error in eval(expr, envir, enclos): 找不到对象'atac_1k'
```

```r
y = atac_1k$TSS.enrichment
```

```
## Error in eval(expr, envir, enclos): 找不到对象'atac_1k'
```

```r
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
```

```
## Error in `[.data.frame`(d, control.list, ): 找不到对象'control.list'
```

```r
#### set a low nucleosome_signal value for 1k stage because it has a great influence on batch effect. 
atac_1k$nucleosome_group <- ifelse(atac_1k$nucleosome_signal > 3 , 'NS > 3' , "NS < 3")
```

```
## Error in ifelse(atac_1k$nucleosome_signal > 3, "NS > 3", "NS < 3"): 找不到对象'atac_1k'
```

```r
FragmentHistogram(atac_1k , group.by = "nucleosome_group")
```

```
## Error in is.data.frame(x): 找不到对象'atac_1k'
```

```r
atac_1k <- RunTFIDF(atac_1k)
```

```
## Error in RunTFIDF(atac_1k): 找不到对象'atac_1k'
```

```r
atac_1k <- FindTopFeatures(atac_1k, min.cutoff = 'q0')
```

```
## Error in FindTopFeatures(atac_1k, min.cutoff = "q0"): 找不到对象'atac_1k'
```

```r
atac_1k <- RunSVD(atac_1k)
```

```
## Error in RunSVD(atac_1k): 找不到对象'atac_1k'
```

```r
atac_1k <- RunUMAP(object = atac_1k, reduction = 'lsi', dims = 2:30)
```

```
## Error in RunUMAP(object = atac_1k, reduction = "lsi", dims = 2:30): 找不到对象'atac_1k'
```

```r
atac_1k <- FindNeighbors(object = atac_1k, reduction = 'lsi', dims = 2:30)
```

```
## Error in FindNeighbors(object = atac_1k, reduction = "lsi", dims = 2:30): 找不到对象'atac_1k'
```

```r
atac_1k <- FindClusters(object = atac_1k, verbose = FALSE, algorithm = 3)
```

```
## Error in FindClusters(object = atac_1k, verbose = FALSE, algorithm = 3): 找不到对象'atac_1k'
```

```r
DimPlot(object = atac_1k, label = TRUE) + NoLegend()
```

```
## Error in is(x, "classRepresentation"): 找不到对象'atac_1k'
```

```r
DimPlot(object = atac_1k, label = TRUE , group.by = "batch")
```

```
## Error in is(x, "classRepresentation"): 找不到对象'atac_1k'
```

```r
DimPlot(object = atac_1k, label = TRUE , group.by = "plate")
```

```
## Error in is(x, "classRepresentation"): 找不到对象'atac_1k'
```

```r
FeaturePlot(atac_1k , features = c("nucleosome_signal" , "TSS.enrichment"))
```

```
## Error in is(x, "classRepresentation"): 找不到对象'atac_1k'
```

```r
######## remove cells that are outliers for these QC metrics.
atac_1k <- subset(atac_1k , 
                    subset = uniq_nuc_frags > 1000 & 
                      TSS.enrichment > 1 & 
                      nucleosome_signal < 3 &
                      mapping_rate > 70)
```

```
## Error in h(simpleError(msg, call)): 在为'subset'函数选择方法时评估'x'参数出了错: 找不到对象'atac_1k'
```

```r
atac_1k <- RunTFIDF(atac_1k)
```

```
## Error in RunTFIDF(atac_1k): 找不到对象'atac_1k'
```

```r
atac_1k <- FindTopFeatures(atac_1k, min.cutoff = 'q0')
```

```
## Error in FindTopFeatures(atac_1k, min.cutoff = "q0"): 找不到对象'atac_1k'
```

```r
atac_1k <- RunSVD(atac_1k)
```

```
## Error in RunSVD(atac_1k): 找不到对象'atac_1k'
```

```r
atac_1k <- RunUMAP(object = atac_1k, reduction = 'lsi', dims = 2:50)
```

```
## Error in RunUMAP(object = atac_1k, reduction = "lsi", dims = 2:50): 找不到对象'atac_1k'
```

```r
DimPlot(object = atac_1k, label = TRUE) + NoLegend()
```

```
## Error in is(x, "classRepresentation"): 找不到对象'atac_1k'
```

```r
DimPlot(object = atac_1k, label = TRUE , group.by = "batch")
```

```
## Error in is(x, "classRepresentation"): 找不到对象'atac_1k'
```

```r
DimPlot(object = atac_1k, label = TRUE , group.by = "plate")
```

```
## Error in is(x, "classRepresentation"): 找不到对象'atac_1k'
```

```r
saveRDS(object = atac_1k , file = "~/zbf_scATAC_n/atac_1k_qc.RDS")
```

```
## Error in saveRDS(object = atac_1k, file = "~/zbf_scATAC_n/atac_1k_qc.RDS"): 找不到对象'atac_1k'
```

This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:


```r
summary(cars)
```

```
##      speed           dist       
##  Min.   : 4.0   Min.   :  2.00  
##  1st Qu.:12.0   1st Qu.: 26.00  
##  Median :15.0   Median : 36.00  
##  Mean   :15.4   Mean   : 42.98  
##  3rd Qu.:19.0   3rd Qu.: 56.00  
##  Max.   :25.0   Max.   :120.00
```

## Including Plots

You can also embed plots, for example:

![plot of chunk pressure](figure/pressure-1.png)

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
