
PolII <- readRDS("~/data/GenomicHeatMap/inst/exdata/PolII.rds")
PolII <- PolII %>% as.data.frame() %>% GenomicRanges::makeGRangesFromDataFrame()
saveRDS(PolII,file = '~/data/GenomicHeatMap/inst/exdata/PolII.rds',compress = 'xz')
rm(PolII)

