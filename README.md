# Zebrafish_scATAC
Scripts and analysis for zebrafish embryo scATAC data in Xu et al.2023, 'The chromatin accessibility dynamics during cell fate specifications in zebrafish early embryogenesis'.

The code provided here is a refinement and supplement to Methods described in our work. Note that this repository is not a container to execute the entire workflow at the paper.

# Overview
The analysis is the downstream of Data preprocessing in paper.

## Data
The outputs of Data preprocessing from each stage were provided and stored in this directory.

## Script
We split the analysis into 7 separate sections and the codes used to reproduce can be found in corresponding R markdown. Files were sorted and numbered in the order described in Methods. To ensure that you would not troubled by 'object not found', we recommend running the code in order, although it is not necessary.

It is worth noting that due to the high degree of repeatability in the following analysis process, we did not fully present the code in order to simplify the space:

**Zf_Scatac_ 1_CreatSignacObject** :  Here we only showed the code used to analyze the 1k stage. The analysis process of samples from other stages is consistent with 1k, so before you run `Zf_ Scatac_ 2_Before IntergrateAllStage`, please refer to the code of `Zf_Scatac_ 1_CreatSignacObject`  to create objects for other stages.

**Zf_scatac_3_NMF** : In the peak annotation section, we only showed the code for module2. When annotating other modules, simply replace 'module2' with the corresponding module.

The version of packages we used can be found in `sessionInfo.md`.

# Contact
If you have any question regarding the code, please open an issue or contact me.

Yunlong Zhang

12131328@mail.sustech.edu.cn
