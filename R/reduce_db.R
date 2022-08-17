
PolII <- readRDS("~/data/GenomicHeatMap/inst/exdata/PolII.rds")
PolII <- PolII %>% as.data.frame() %>% GenomicRanges::makeGRangesFromDataFrame()
saveRDS(PolII,file = '~/data/GenomicHeatMap/inst/exdata/PolII.rds',compress = 'xz')
rm(PolII)

H3K79me2 <- readRDS("~/data/GenomicHeatMap/inst/exdata/H3K79me2.rds")
H3K79me2 <- H3K79me2 %>% as.data.frame() %>% GenomicRanges::makeGRangesFromDataFrame()
saveRDS(H3K79me2,file = '~/data/GenomicHeatMap/inst/exdata/H3K79me2.rds',compress = 'xz')
