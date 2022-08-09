
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GenomicHeatMap

<!-- badges: start -->
<!-- badges: end -->

GenomicHeatMap is used to generate genomics and epigeniomic heatmaps
from insertion sites

## Installation

You can install the development version of GenomicHeatMap from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Adrian-Cantu/GenomicHeatMap")
```

## Example

You can get an epigenomic heatmap plot,

``` r
library(GenomicHeatMap)
#intsite_to_heatmap_df(intSites) %>% epi_annotate_df() -> hh
#hh %>% make_roc() %>% make_heatmap(title='nice heatmap')
```
