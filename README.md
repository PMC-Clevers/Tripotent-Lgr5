
#  Tripotent Lgr5 stem cells in the posterior tongue generate non-taste lingual, taste bud and salivary gland cells

This repository contains R scripts used to generate **Figures 1,2,4,6,7** from the manuscript. It focuses on dissecting the transcriptional diversity and developmental trajectories of taste bud and salivary gland lineages using single-cell RNA sequencing (scRNA-seq).

Data and code availability 

Processed bulk and single-cell RNA sequencing data are deposited to the GEO repository(GSE274014(Foxe1 knock-out RNA-seq), GSE274015 (scRNA-seq clonal and tissue)  and are available in a database with controlled access. Code will be made available upon publication on https://github.com/PMC-Clevers/Tripotent-Lgr5/.

---

## 📁 Project Overview

| Figure | Description |
|--------|-------------|
| **Figure 1** | Single cell atlas projects the complete vivid tapestry of the posterior tongue |
| **Figure 2** | Characterization of the taste buds and minor salivary glands on the posterior tongue |
| **Figure 4** | Lgr5 stem cell is capable of forming non-taste lingual epithelium, salivary gland and taste bud lineages|
| **Figure 6** | Transcriptional characterization of cell fate decisions|
| **Figure 7** | Foxe1-KO organoids project an ablation of differentiation queues|
| **Supplementary Figures** | Marker gene expression (DotPlots), GO enrichment, and von Ebner gland analysis |

NOTE: scripts for easy featurplots showing relative expression are not projected in the code as they are simple plots
---

## 📦 Datasets Used

- `cvp`, `fop`: Raw count matrices from CVP and FoP
- `annocvpfop`: Annotated merged Seurat object for taste tissues
- `postnatalSMG`: External dataset for acinar subtype annotation
- `vonEbner`: Subset of annocvpfop corresponding to von Ebner gland
- `taste`: Subset of annocvpfop corresponding to von Ebner gland
- Clone datasets: `Clone1`, `Clone3`, `Clone7`

---

## 🧪 Required Packages

All figures rely on the following R packages:

```r
library(Signac)
library(Seurat)
library(SeuratData)
library(scater)
library(SeuratWrappers)
library(monocle3)
library(monocle)
library(Matrix)
library(ggplot2)
library(patchwork)
library(GeneSwitches)
library(SingleCellExperiment)
library(dplyr)
library(viridis)
library(RColorBrewer)
library(pheatmap)
library(circlize)
library(scCustomize)
library(clusterProfiler)
library(org.Mm.eg.db)
library(EnhancedVolcano)
library(ggplot2)
library(scales) # needed for oob parameter
library(viridis)
library(reshape2)
library(fgsea)
library(DESeq2)
library(dplyr)
library(pheatmap)


sessionInfo()
R version 4.4.0 (2024-04-24)
Platform: x86_64-apple-darwin20
Running under: macOS 15.3.2

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.4-x86_64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: Europe/Amsterdam
tzcode source: internal

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] scCustomize_3.0.1           circlize_0.4.16             RColorBrewer_1.1-3         
 [4] GeneSwitches_0.1.0          patchwork_1.3.0             Matrix_1.7-0               
 [7] monocle3_1.3.7              SeuratWrappers_0.3.5        scater_1.32.1              
[10] scuttle_1.14.0              SingleCellExperiment_1.26.0 SeuratData_0.2.2.9001      
[13] Seurat_5.2.1                SeuratObject_5.0.2          sp_2.2-0                   
[16] Signac_1.14.0               pheatmap_1.0.12             dplyr_1.1.4                
[19] DESeq2_1.44.0               SummarizedExperiment_1.34.0 Biobase_2.64.0             
[22] MatrixGenerics_1.16.0       matrixStats_1.5.0           GenomicRanges_1.56.2       
[25] GenomeInfoDb_1.40.1         IRanges_2.38.1              S4Vectors_0.42.1           
[28] BiocGenerics_0.50.0         fgsea_1.30.0                reshape2_1.4.4             
[31] viridis_0.6.5               viridisLite_0.4.2           scales_1.3.0               
[34] EnhancedVolcano_1.22.0      ggrepel_0.9.6               ggplot2_3.5.1              

loaded via a namespace (and not attached):
  [1] RcppAnnoy_0.0.22          splines_4.4.0             later_1.4.1              
  [4] bitops_1.0-9              R.oo_1.27.0               tibble_3.2.1             
  [7] polyclip_1.10-7           janitor_2.2.1.9000        fastDummies_1.7.5        
 [10] lifecycle_1.0.4           Rdpack_2.6.2              globals_0.16.3           
 [13] lattice_0.22-6            MASS_7.3-60.2             magrittr_2.0.3           
 [16] plotly_4.10.4             remotes_2.5.0             httpuv_1.6.15            
 [19] sctransform_0.4.1         spam_2.11-1               spatstat.sparse_3.1-0    
 [22] reticulate_1.41.0.1       cowplot_1.1.3             pbapply_1.7-2            
 [25] minqa_1.2.8               lubridate_1.9.4           abind_1.4-8              
 [28] zlibbioc_1.50.0           Rtsne_0.17                mixtools_2.0.0.1         
 [31] R.utils_2.13.0            purrr_1.0.4               RCurl_1.98-1.16          
 [34] rappdirs_0.3.3            GenomeInfoDbData_1.2.12   irlba_2.3.5.1            
 [37] listenv_0.9.1             spatstat.utils_3.1-2      goftest_1.2-3            
 [40] RSpectra_0.16-2           bigmemory_4.6.4           spatstat.random_3.3-2    
 [43] fitdistrplus_1.2-2        parallelly_1.42.0         DelayedMatrixStats_1.26.0
 [46] codetools_0.2-20          DelayedArray_0.30.1       RcppRoll_0.3.1           
 [49] shape_1.4.6.1             tidyselect_1.2.1          UCSC.utils_1.0.0         
 [52] farver_2.1.2              lme4_1.1-36               ScaledMatrix_1.12.0      
 [55] spatstat.explore_3.3-4    jsonlite_1.9.1            BiocNeighbors_1.22.0     
 [58] progressr_0.15.1          ggridges_0.5.6            survival_3.6-4           
 [61] segmented_2.1-4           tools_4.4.0               ica_1.0-3                
 [64] Rcpp_1.0.14               glue_1.8.0                gridExtra_2.3            
 [67] SparseArray_1.4.8         withr_3.0.2               BiocManager_1.30.25      
 [70] fastmap_1.2.0             boot_1.3-30               rsvd_1.0.5               
 [73] digest_0.6.37             timechange_0.3.0          R6_2.6.1                 
 [76] mime_0.12                 ggprism_1.0.5             colorspace_2.1-1         
 [79] scattermore_1.2           tensor_1.5                spatstat.data_3.1-4      
 [82] R.methodsS3_1.8.2         tidyr_1.3.1               generics_0.1.3           
 [85] data.table_1.17.0         httr_1.4.7                htmlwidgets_1.6.4        
 [88] S4Arrays_1.4.1            uwot_0.2.3                pkgconfig_2.0.3          
 [91] gtable_0.3.6              fastglm_0.0.3             lmtest_0.9-40            
 [94] XVector_0.44.0            htmltools_0.5.8.1         dotCall64_1.2            
 [97] png_0.1-8                 snakecase_0.11.1          spatstat.univar_3.1-2    
[100] reformulas_0.4.0          bigmemory.sri_0.1.8       rstudioapi_0.17.1        
[103] uuid_1.2-1                nlme_3.1-164              nloptr_2.2.0             
[106] GlobalOptions_0.1.2       zoo_1.8-13                stringr_1.5.1            
[109] KernSmooth_2.23-24        vipor_0.4.7               parallel_4.4.0           
[112] miniUI_0.1.1.1            ggrastr_1.0.2             pillar_1.10.1            
[115] grid_4.4.0                vctrs_0.6.5               RANN_2.6.2               
[118] promises_1.3.2            BiocSingular_1.20.0       beachmat_2.20.0          
[121] xtable_1.8-4              cluster_2.1.6             paletteer_1.6.0          
[124] beeswarm_0.4.0            cli_3.6.4                 locfit_1.5-9.12          
[127] compiler_4.4.0            Rsamtools_2.20.0          rlang_1.1.5              
[130] crayon_1.5.3              future.apply_1.11.3       rematch2_2.1.2           
[133] forcats_1.0.0             ggbeeswarm_0.7.2          plyr_1.8.9               
[136] stringi_1.8.4             deldir_2.0-4              BiocParallel_1.38.0      
[139] munsell_0.5.1             Biostrings_2.72.1         lazyeval_0.2.2           
[142] spatstat.geom_3.3-5       RcppHNSW_0.6.0            sparseMatrixStats_1.16.0 
[145] future_1.34.0             shiny_1.10.0              kernlab_0.9-33           
[148] rbibutils_2.3             ROCR_1.0-11               igraph_2.1.4             
[151] fastmatch_1.1-6     
