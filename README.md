# Zebrafish_scATAC
Scripts and analysis for zebrafish embryo scATAC data in Xu et al.2023, [The chromatin accessibility dynamics during cell fate specifications in zebrafish early embryogenesis](https://www.biorxiv.org/content/10.1101/2023.10.13.562312v1).

The code provided here is a refinement and supplement to Methods described in our work. Note that this repository is not a container to execute the entire workflow at the paper.

# Overview
The analysis is the downstream of Data preprocessing in paper.

## Data
Meta data after data preprocess and filter was saved in `metadata.csv`, including batches and cell types.

## Script
We split the analysis into 7 separate sections and the codes used to reproduce can be found in corresponding R markdown. Files were sorted and numbered in the order described in Methods. To ensure that you would not troubled by 'object not found', we recommend running the code in order, although it is not necessary. 

It is worth noting that due to the high degree of code duplication during the following analysis process, in order to avoid verbosity, we did not fully present the code:

**Creat Signac Object** :  We only showed the codes used to create object for data from **3.0 hpf** in `Zf_scatac_ 1_CreatSignacObject.Rmd`. The process to create object for data from other stages is consistent with it, so before run `Zf_ scatac_ 2_ IntergrateAllStage.Rmd`, please refer to the codes of `Zf_scatac_ 1_CreatSignacObject.Rmd` to create objects for other stages.

**Peaks Annotation** : In the peak annotation section of `Zf_scatac_3_NMF.Rmd`, we only showed the codes for module2. When annotating other modules, simply replace 'module2' with the corresponding module.

The version of packages we used can be found in `sessionInfo.md`.

# Contact
If you have any question regarding the code, please open an issue or contact me.

Yunlong Zhang

12131328@mail.sustech.edu.cn
