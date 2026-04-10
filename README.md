# README

This repository contains the code and data for:

**Correcting STKDE Using SPD: A Case Study from the Middle Neolithic of Northeastern Iberia**

## Paper Authors

Paper Authors Biel Soriano Elias (a), Andreu Monforte-Barberán (a), Anna Bach Gómez (a) & Miquel Molist (a)

a Autonomous University of Barcelona, Prehistory Departament, SAPPO-GRAMPO

## Structure of the repository

The repository folder is structured as follows:

- **README.md**: This file (repository overview).  
- **Code/**: Contains all R scripts for the paper.  
  - **Main **: Code (R scripts) to conduct all the analysis in the present paper
- **Data/**: Raw data used in the paper, all in the CRS ETRS89 / UTM 31 N.  
  - **Rasters/**: All the necessary rasters for the paper, including DEM, land-use, RBIAS and combined rasters.
    - **MDE_cat_and_100**: DEM (.tiff) cropped and reprojected from the GLO-30 Copernicus
    (https://ec.europa.eu/eurostat/web/gisco/geodata/digital-elevation-model/copernicus#Elevation, last accessed on 9/1/2026)
  - **RAW_burials/**: Data of burials
    - **c14_raw_burials.xlsx**: Database of C14 dated burials (.xlsx)
    - **no_c14_raw_burials.xlsx**: Database of undated burials (.xlsx)
    
## Computational Environment

All analyses and code development were conducted on:

- Windows 10 (64-bit): HP EliteBook 640 14inch G9 Notebook, 12th Gen Intel(R) Core (TM) i5, 16GB RAM  

### Information about the R Session

R version 4.5.3 (2026-03-11)

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] reshape2_1.4.4         sqldf_0.4-12           RSQLite_2.4.6          gsubfn_0.7             proto_1.0.0            BBmisc_1.13.1         
 [7] Bchron_4.7.8           data.table_1.17.8      gdata_3.0.1            scales_1.4.0           Metrics_0.1.4          vegan_2.7-3           
[13] permute_0.9-10         mclust_6.1.2           corrplot_0.95          magick_2.9.1           rcarbon_1.5.2          viridis_0.6.5         
[19] viridisLite_0.4.2      stars_0.6-8            abind_1.4-8            spatstat_3.5-1         spatstat.linnet_3.4-1  spatstat.model_3.6-1  
[25] rpart_4.1.24           spatstat.explore_3.7-0 nlme_3.1-168           spatstat.random_3.4-4  spatstat.geom_3.7-0    spatstat.univar_3.1-6 
[31] spatstat.data_3.1-9    ggplot2_4.0.2          doParallel_1.0.17      iterators_1.0.14       foreach_1.5.2          tidyr_1.3.2           
[37] dplyr_1.2.0            sf_1.0-24              terra_1.8-93           readxl_1.4.5          

loaded via a namespace (and not attached):
 [1] DBI_1.2.3             deldir_2.0-4          gridExtra_2.3         tcltk_4.5.3           rlang_1.1.7           magrittr_2.0.3       
 [7] ggridges_0.5.7        e1071_1.7-16          compiler_4.5.3        mgcv_1.9-4            vctrs_0.7.1           stringr_1.5.1        
[13] pkgconfig_2.0.3       fastmap_1.2.0         backports_1.5.0       rmarkdown_2.29        purrr_1.0.4           bit_4.6.0            
[19] xfun_0.52             cachem_1.1.0          goftest_1.2-3         blob_1.3.0            spatstat.utils_3.2-1  tweenr_2.0.3         
[25] cluster_2.1.8.2       R6_2.6.1              stringi_1.8.7         RColorBrewer_1.1-3    cellranger_1.1.0      Rcpp_1.0.14          
[31] knitr_1.50            tensor_1.5            snow_0.4-4            Matrix_1.7-4          splines_4.5.3         tidyselect_1.2.1     
[37] rstudioapi_0.17.1     yaml_2.3.10           codetools_0.2-20      plyr_1.8.9            lattice_0.22-9        tibble_3.2.1         
[43] withr_3.0.2           S7_0.2.1              evaluate_1.0.5        units_0.8-7           proxy_0.4-27          polyclip_1.10-7      
[49] pillar_1.11.0         KernSmooth_2.23-26    checkmate_2.3.4       generics_0.1.4        chron_2.3-62          gtools_3.9.5         
[55] class_7.3-23          glue_1.8.0            tools_4.5.3           grid_4.5.3            ggforce_0.5.0         cli_3.6.5            
[61] spatstat.sparse_3.1-0 doSNOW_1.0.20         gtable_0.3.6          digest_0.6.37         classInt_0.4-11       farver_2.1.2         
[67] memoise_2.0.1         htmltools_0.5.8.1     lifecycle_1.0.5       bit64_4.6.0-1         MASS_7.3-65          