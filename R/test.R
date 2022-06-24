# library
library(tidyverse)
library(hiAnnotator)
library(purrr)
library(hotROCs)
library(intSiteRetriever)
get_N_MRCs()

# kk <- get_reference_genome('hg38')

genome_sequence <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
genome_sequence@user_seqnames <- genome_sequence@user_seqnames[genome_sequence@user_seqnames %in% paste0("chr", c(1:22, "X", "Y", "M"))]
genome_sequence@seqinfo <- genome_sequence@seqinfo[paste0("chr", c(1:22, "X", "Y", "M"))]

intSites <- readRDS('intSites.rds')

intSites_coor <- intSites %>%
  as.data.frame() %>%
  select(c(seqnames,start,end,strand,GTSP)) %>%
  mutate(type='insertion')



# df_to_randomize <- intSites %>%
#   as.data.frame() %>%
#   select(c(seqnames,strand,start)) %>%
#   dplyr::rename(position=start) %>%
#   dplyr::rename(chr=seqnames) %>%
#   mutate(siteID=paste0('site',row_number())) %>%
#   relocate(siteID)

df_to_randomize <- tibble(.rows = length(intSites)) %>%
  mutate(siteID=row_number()) %>%
  mutate(gender='m')



set.seed(as.numeric(Sys.time()))
ttt<-get_N_MRCs(df_to_randomize,genome_sequence)

tttt <- ttt %>%
  dplyr::rename(start=position) %>%
  mutate(end=start) %>%
  mutate(siteID=NULL) %>%
  relocate(strand,.after = last_col()) %>%
  mutate(GTSP=rep(intSites_coor$GTSP,3)) %>%
  mutate(type='match') %>%
  dplyr::rename(seqnames=chr)

to_get_features <- rbind(intSites_coor,tttt)


epi_files <- list.files('epigenetic_features_d')
names(epi_files) <- epi_files %>% str_remove(., ".rds")


#all_names <- unique(intSites$GTSP)


all_epi <- lapply(epi_files,function(x){readRDS(file.path('epigenetic_features_d', x))})
#kk_test <- getFeatureCounts(intSites, all_epi[[1]], 'test')


  kk <- imap(all_epi, function(x,name){
    to_get_features <<- getFeatureCounts(to_get_features, x, name)
    return(1)
  })

#head(to_get_features)

to_roc_df <-as.data.frame(to_get_features)

# get_annotation_columns <- function(sites) {
#   granges_column_names <- c("seqnames", "start", "end", "width", "strand")
#   int_site_column_names <- c("siteID", "sampleName", "chr", "strand", "position")
#   required_columns <- unique(c(
#     granges_column_names, int_site_column_names, "type"))
#   stopifnot(all(required_columns %in% names(sites)))
#   setdiff(names(sites), required_columns)
# }
# setdiff(colnames(jointsites), required_columns)
# required_columns
# required_columns %in% colnames(jointsites)
# get_annotation_columns(jointsites)

roc.res <- ROC.ORC(
  response = to_roc_df$type,
  variables = to_roc_df %>% select(-c(seqnames,start,end,width,strand,GTSP,type)),
  origin=to_roc_df$GTSP)

roc.res$ROC %>%
  as.data.frame() %>%
  rownames_to_column(var = "feature") %>%
  pivot_longer(!feature,values_to='val', names_to='sample') %>%
  ggplot( aes(sample,feature, fill= val)) +
  geom_tile() +
  scale_fill_gradientn(colours=c('blue','white','red'),
                       na.value = "transparent",
                       breaks=c(0,0.5,1),
                       labels=c(0,0.5,1),
                       limits=c(0,1))



