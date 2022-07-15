# library
library(tidyverse)
library(hiAnnotator)
library(purrr)
library(hotROCs)
library(intSiteRetriever)
#get_N_MRCs()
library(GenomicRanges)
#jurkat <- read_tsv('wgEncodeRegDnaseUwJurkatPeak.tsv')
#jurkat_gr <- jurkat %>% mutate(`#bin`=NULL) %>% makeGRangesFromDataFrame()
#saveRDS(jurkat_gr,file=file.path('epigenetic_features_d','DnaseUwJurkat.rds'))
# kk <- get_reference_genome('hg38')




genome_sequence <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
genome_sequence@user_seqnames <- genome_sequence@user_seqnames[genome_sequence@user_seqnames %in% paste0("chr", c(1:22, "X", "Y", "M"))]
genome_sequence@seqinfo <- genome_sequence@seqinfo[paste0("chr", c(1:22, "X", "Y", "M"))]


# DNASE1 <- read_table('wgEncodeRegDnaseClustered.txt',col_names = c('bin','chrom','start','end','name','score',
#                                                                    'sourceCounts','sourceIds','sourceScores'))
# DNASE1_gr <- DNASE1 %>% filter(chrom %in% paste0("chr", c(1:22, "X", "Y", "M"))) %>%
#   select(c(chrom,start,end)) %>%
#   makeGRangesFromDataFrame(seqinfo=genome_sequence@seqinfo)
# saveRDS(DNASE1_gr,file=file.path('epigenetic_features_d','Dnase.rds'))

# CPGis <- read_table('cpgIslandExt.txt',
#                     col_names = c('bin','chrom','start','end','v1','v2','v3','v4','v5','v6','v7','v8'))
# CPGis_gr <- CPGis %>% filter(chrom %in% paste0("chr", c(1:22, "X", "Y", "M"))) %>%
#   select(c(chrom,start,end)) %>%
#   makeGRangesFromDataFrame(seqinfo=genome_sequence@seqinfo)
#   saveRDS(CPGis_gr,file=file.path('epigenetic_features_d','CPGis.rds'))

library(RMySQL)
library(gt23)

histone_roc_features <- c('H3K4me1',
                          'H3K4me2',
                          'H3K4me3',
                          'H3K9me1',
                          'H3K9me2',
                          'H3K9me3',
                          'H3K27me1',
                          'H3K27me2',
                          'H3K27me3',
                          'H3K36me1',
                          'H3K36me3',
                          'H3K79me1',
                          'H3K79me2',
                          'H3K79me3',
                          'H3R2me1',
                          'H3R2me2',
                          'H4K20me1',
                          'H4K20me3',
                          'H3K4ac',
                          'H3K9ac',
                          'H3K14ac',
                          'H3K18ac',
                          'H3K23ac',
                          'H3K27ac',
                          'H3K36ac',
                          'H4K5ac',
                          'H4K8ac',
                          'H4K12ac',
                          'H4K16ac',
                          'H4K91ac',
                          'H2AK5ac',
                          'H2AK9ac',
                          'H2BK5ac',
                          'H2BK12ac',
                          'H2BK20ac',
                          'H2BK120ac',
                          'H2BK5me1',
                          'H2A.Z'
)

intSites <- readRDS("20220714_ViiV_insites_plus_sampleinfo.rds")

### need to chose the right grouping variable, maybe patient or GTSP or some other
intSites_coor <- intSites %>%
  as.data.frame() %>%
  filter(to_heatmap) %>%
  select(c(seqnames,start,end,strand,patient,GTSP,Drug,concentration_nM,concentration_txt,Replicate)) %>%
  mutate(type='insertion')



# df_to_randomize <- intSites %>%
#   as.data.frame() %>%
#   select(c(seqnames,strand,start)) %>%
#   dplyr::rename(position=start) %>%
#   dplyr::rename(chr=seqnames) %>%
#   mutate(siteID=paste0('site',row_number())) %>%
#   relocate(siteID)




epi_files <- list.files('epigenetic_features_d')
names(epi_files) <- epi_files %>% str_remove(., ".rds")


#all_names <- unique(intSites$GTSP)

### need to fix this no it is not all loaded in ram
#all_epi <- lapply(epi_files,function(x){readRDS(file.path('epigenetic_features_d', x))})
#kk_test <- getFeatureCounts(intSites, all_epi[[1]], 'test')


#to_roc_df %>%  filter(if_any(everything(), is.na)) %>% replace(is.na(.), 0)

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

if(! file.exists('roc_res.rds')){
  df_to_randomize <- tibble(.rows = nrow(intSites_coor)) %>%
    mutate(siteID=row_number()) %>%
    mutate(gender='m')



  set.seed(as.numeric(Sys.time()))
  ttt<-get_N_MRCs(df_to_randomize,genome_sequence)

  #set how the groups are defined by changeing patient
  group_vector<-paste0(intSites_coor$Drug,'_',intSites_coor$concentration_txt)
  intSites_coor$patient <- group_vector


  tttt <- ttt %>%
    dplyr::rename(start=position) %>%
    mutate(end=start) %>%
    mutate(siteID=NULL) %>%
    relocate(strand,.after = last_col()) %>%
    mutate(patient=rep(group_vector,3)) %>%
    mutate(type='match') %>%
    dplyr::rename(seqnames=chr)

  to_get_features <- rbind(intSites_coor %>% select(colnames(tttt)) ,tttt)

  kk <- imap(epi_files, function(x,name){
    xx <- readRDS(file.path('epigenetic_features_d', x))
    to_get_features <<- getFeatureCounts(to_get_features, xx, name)
    rm(xx)
    gc()
    return(1)
  })

  #head(to_get_features)

  to_roc_df <-as.data.frame(to_get_features) %>% replace(is.na(.), 0)


roc.res <- ROC.ORC(
  response = to_roc_df$type,
  variables = to_roc_df %>% select(-c(seqnames,start,end,width,strand,patient,type)),
  origin=to_roc_df$patient)

saveRDS(roc.res,file = 'roc_res.rds')
} else {
  roc.res <- readRDS(file='roc_res.rds')
}


sort_features <- function(...){

  translate_window <- setNames(c(1000,1000000), c("Kb", "Mb"))
  xx <- tibble(.rows = length(unique(as.vector(...)))) %>%
    mutate(site=unique(as.vector(...))) %>%
    separate(site,into=c('fname','window'),sep = '\\.',remove = FALSE) %>%
    separate(window,into=c('size','todel'),sep = '\\D+$',remove = FALSE) %>%
    separate(window,into=c('todel','size_word'),sep = '^\\d+',remove = FALSE) %>%
    mutate(size_mult=translate_window[size_word]) %>%
    mutate(size_sort=as.numeric(size)*size_mult) %>%
    mutate(todel=NULL) %>%
    arrange(fname,rev(size_sort)) %>%
    pull(var = site)

  toret <- factor(as.vector(...),xx)
  return(toret)
}

roc_df <- roc.res$ROC %>%
  as.data.frame() %>%
  rownames_to_column(var = "feature") %>%
  pivot_longer(!feature,values_to='val', names_to='sample') %>%
  mutate(feature=sort_features(feature)) %>%
  separate(feature,into = c('feature_name','feature_concentration'),sep = '\\.',remove = FALSE)

stat_cuts <- c(0, 0.001, 0.01, 0.05, 1)
roc_pval <- roc.res$pvalues$np %>%
  as.data.frame() %>%
  rownames_to_column(var = "feature") %>%
  pivot_longer(!feature,values_to='pval', names_to='sample') %>%
  mutate(feature=sort_features(feature)) %>%
  mutate(pval_txt=cut(pval, stat_cuts, labels = c("***", " **", " * ", "   "),include.lowest = TRUE))

his_roc_df <- left_join(roc_df %>% filter(feature_name %in% histone_roc_features) %>%
                          mutate(jj=paste0(feature,'_',sample)),
          roc_pval%>% mutate(jj=paste0(feature,'_',sample)) %>% select(c(pval,pval_txt,jj)),
          by='jj') %>%
  mutate(jj=NULL)


his_roc_df %>%
  ggplot( aes(y=feature,x=sample, fill= val)) +
  geom_tile() +
  geom_text(aes(label = pval_txt), color = "black", size = 3, nudge_y = -0.15)+
  theme_classic() +
  #theme_minimal()+
  scale_y_discrete(labels=his_roc_df$feature_concentration,breaks=his_roc_df$feature)+
  #geom_text(position = position_dodge(width = 1), aes(x=0, y=feature, label=feature_name)) +
  #facet_wrap(~feature_name, scales="free")+

  labs(fill="ROC area",title="Histone heatmap") +
  scale_fill_gradientn(colours=c('blue','grey90','red'),
                       na.value = "transparent",
                       breaks=c(0,0.5,1),
                       labels=c(0,0.5,1),
                       limits=c(0,1)) +
  theme(axis.text.x = element_text(angle = 30,vjust = 1, hjust=1),
        axis.text.y.left = element_text(size=7),
        axis.title.x=element_blank(),
        panel.spacing.y = unit(-0.15, "line"),
        strip.placement='outside',
        panel.border = element_blank(),
        panel.background= element_blank(),
        strip.background = element_blank(),
        strip.text.y.left = element_text(angle=0,size=10),
        plot.background = element_blank(),
        # legend.position="top",
        axis.ticks.length.y.left=unit(0.1,'line'),
        axis.title.y=element_blank()) +
  facet_grid(rows=vars(his_roc_df$feature_name),scales = "free_y",switch = 'y')+
#  coord_fixed(ratio=2/1) +
  scale_x_discrete(expand = c(0,0)) -> pp



library(ggplot2)
library(grid)
#ggsave('viiv_epi_genmap_test_20220712.pdf',width = 8.5,height = 11,units = 'in')

#data = data.frame(feature_name='H4K16ac')

#p1 = pp + annotation_custom2(linesGrob(), xmin = 2.6, xmax = 2.75, ymin = 100, ymax = 100, data = his_roc_df)

#tag_facet3(pp,x1=11,x2=11,y1=-5,y2=-10)

pdf(file='histone_heatmap.pdf',width = 8.5,height = 11)
#grid::grid.newpage()
#vp <- viewport(width = unit(7.5, "inches"), height = unit(10, "inches"))
#pushViewport(vp)
grid::grid.draw(pp)
y1<-0.07
dy<-0.0184
ddy<-0.006
dx <- 0.01
grid.lines(x=c(0.11,0.11),y=c(y1,y1+dy),gp=grid::gpar(lwd=1))
grid.lines(x=c(0.11,0.11+dx),y=c(y1,y1),gp=grid::gpar(lwd=1))
grid.lines(x=c(0.11,0.11+dx),y=c(y1+dy,y1+dy),gp=grid::gpar(lwd=1))
for (val in 1:36) {
  y1<-y1+ddy+dy
  grid.lines(x=c(0.11,0.11),y=c(y1,y1+dy),gp=grid::gpar(lwd=1))
  grid.lines(x=c(0.11,0.11+dx),y=c(y1,y1),gp=grid::gpar(lwd=1))
  grid.lines(x=c(0.11,0.11+dx),y=c(y1+dy,y1+dy),gp=grid::gpar(lwd=1))
}
dev.off()

#ggsave('test.pdf',width = 8.5,height = 11,units = 'in')
# # Code to override clipping
# gt <- ggplotGrob(p1)
# gt$layout[grepl("panel", gt$layout$name), ]$clip <- "off"
#
# # Draw the plot
# grid.newpage()
# grid.draw(gt)



  #separate(feature,into = c('feature_name','feature_concentration'),sep = '\\.',remove = FALSE)



#roc.res$ROC
