# library
library(tidyverse)
library(hiAnnotator)
library(purrr)
library(hotROCs)
library(intSiteRetriever)
#get_N_MRCs()
library(GenomicRanges)
library(rtracklayer)

#
genome_sequence <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
genome_sequence@user_seqnames <- genome_sequence@user_seqnames[genome_sequence@user_seqnames %in% paste0("chr", c(1:22, "X", "Y", "M"))]
genome_sequence@seqinfo <- genome_sequence@seqinfo[paste0("chr", c(1:22, "X", "Y", "M"))]
#

intSites <- readRDS("20220714_ViiV_insites_plus_sampleinfo.rds")

### need to chose the right grouping variable, maybe patient or GTSP or some other
intSites_coor <- intSites %>%
  as.data.frame() %>%
  filter(to_heatmap) %>%
  select(c(seqnames,start,end,strand,patient,GTSP,Drug,concentration_nM,concentration_txt,Replicate)) %>%
  mutate(type='insertion')


### start roc code
if(! file.exists('roc_res_gen.rds')){
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


#################################### start getting genomic ---------------------

## windows
window_size <- c("1k"=1e4, "10k"=1e4, "1M"=1e6)

#
refGenes <- readRDS(file="hg38.refSeq.rds")
refGenes <- refGenes[seqnames(refGenes) %in% paste0("chr", c(1:22, "X", "Y", "M"))]

CpG_islands <- readRDS(file=file.path('epigenetic_features_d','CPGis.rds'))

DNaseI <- readRDS(file=file.path('epigenetic_features_d','Dnase.rds'))


#
from_counts_to_density <- function(sites, column_prefix, window_size) {
  metadata <- mcols(sites)
  sapply(seq(window_size), function(i) {
    val <- window_size[i]
    name <- names(window_size)[i]
    column_name <- paste0(column_prefix, ".", name)
    metadata[[column_name]] <<- metadata[[column_name]]/val
  })
  mcols(sites) <- metadata
  sites
}


hg38_seqinfo <- genome_sequence@seqinfo



  c_sample <- to_get_features

  gen_final_data <- c_sample %>% #as.data.frame() %>% slice_head(n=12000) %>% slice_tail(n=6000) %>%  # sample_n(40) %>%   #filter(start==142938751) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=TRUE,seqinfo = hg38_seqinfo) %>%
    hiAnnotator::getFeatureCounts(refGenes, "refSeq_counts", width = window_size) %>%
    GCcontent::getGCpercentage("GC", window_size, genome_sequence) %>%
    hiAnnotator::getFeatureCounts(CpG_islands, "CpG_counts", width = window_size) %>%
    hiAnnotator::getFeatureCounts(CpG_islands, "CpG_density", width = window_size) %>%
    from_counts_to_density("CpG_density", window_size) %>%
#    hiAnnotator::getFeatureCounts(DNaseI, "DNaseI_count", width = window_size) %>%
    as.data.frame()

# debug area -----------

# gen_final_data %>% filter(if_any(everything(), is.na)) %>% head()
# testhead <- gen_final_data %>% filter(!is.na(GC.1k)) %>% filter(is.na(GC.100)) %>%
#   GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=TRUE,seqinfo = hg38_seqinfo)
#


#Biostrings::getSeq(genome_sequence,testhead+500)

#gen_final_data %>% filter(if_any(everything(), is.na)) %>% pull(var = type) %>% unique()
#gen_final_data %>% filter(if_any(everything(), is.na)) %>%  head()

# DNAse1 and cpg island where not downloading fully, fixed getting file from UCSC server
# DNAse1 and cpg island are not annotated in ChrmM, just affect 1 site, ignorign for the moment
# CG content retrun NA if there are only N, affect only random samples

###########

to_roc_df <-as.data.frame(gen_final_data) %>% replace(is.na(.), 0)
roc.res <- ROC.ORC(
  response = to_roc_df$type,
  variables = to_roc_df %>% select(-c(seqnames,start,end,width,strand,patient,type)),
  origin=to_roc_df$patient)

saveRDS(roc.res,file = 'roc_res_gen.rds')
} else {
  roc.res <- readRDS(file='roc_res_gen.rds')
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


roc_df_tmp <- roc.res$ROC %>%
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

roc_df_tmp2 <- left_join(roc_df_tmp %>%  mutate(jj=paste0(feature,'_',sample)),
                        roc_pval%>% mutate(jj=paste0(feature,'_',sample)) %>% select(c(pval,pval_txt,jj)),
                        by='jj') %>%
  mutate(jj=NULL)

#### vehicle star

mk.stars <- function(pvmat) {
  x <- array("", dim(pvmat))
  x[] <- as.character(cut(pvmat, c(0, 0.001, 0.01, 0.05,
                                   1), c("***", "**", "*", ""), include.lowest = TRUE))
  x
}

roc.rows <- nrow(roc.res$ROC)
roc.cols <- ncol(roc.res$ROC)

opvals <- roc.res$pvalues$op
vpvals <- roc.res$pvalues$vp
nullpvals <- roc.res$pvalues$np
opstars <- mk.stars(opvals)
vpstars <- mk.stars(vpvals)
nullpstars <- mk.stars(nullpvals)
matchCol <- function(x, indx) x[2:1, ][x == indx]
isCol <- function(x, indx) colSums(x == indx) == 1
omasks <- lapply(1:roc.cols, function(x) {
  res <- array("", dim(roc.res$ROC))
  res[, x] <- "--"
  wc <- attr(opvals, "whichCol")
  res[, matchCol(wc, x)] <- opstars[, isCol(wc, x)]
  res
})

vehicle_stars <- omasks[[7]] %>%
  `colnames<-`(colnames(roc.res$ROC)) %>%
  `rownames<-`(rownames(roc.res$ROC)) %>%
  as.data.frame() %>%
  rownames_to_column(var = "feature") %>%
  pivot_longer(!feature,values_to='v_star', names_to='sample')


roc_df <- left_join(roc_df_tmp2 %>%  mutate(jj=paste0(feature,'_',sample)),
                    vehicle_stars%>% mutate(jj=paste0(feature,'_',sample)) %>% select(c(v_star,jj)),
                    by='jj') %>% mutate(jj=NULL)


####


roc_df %>%
  ggplot( aes(y=feature,x=sample, fill= val)) +
  geom_tile() +
  #geom_text(aes(label = pval_txt), color = "black", size = 3, nudge_y = -0.15)+
  geom_text(aes(label = v_star), color = "black", size = 5, nudge_y = -0.15)+
  theme_classic() +
  scale_y_discrete(labels=roc_df$feature_concentration,breaks=roc_df$feature,expand = c(0,0))+
  labs(fill="ROC area",title="Genomic heatmap") +
  scale_fill_gradientn(colours=c('purple2','grey95','yellow'),
                       na.value = "transparent",
                       breaks=c(0,0.5,1),
                       labels=c(0,0.5,1),
                       limits=c(0,1)) +
  theme(axis.text.x = element_text(angle = 30,vjust = 1, hjust=1),
        axis.text.y.left = element_text(size=11),
        axis.title.x=element_blank(),
        panel.spacing.y = unit(-0.1, "line"),
        strip.placement='outside',
        panel.border = element_blank(),
        panel.background= element_blank(),
        strip.background = element_blank(),
        strip.text.y.left = element_text(angle=0,size=12),
        plot.background = element_blank(),
        axis.ticks.length.y.left=unit(0.1,'line'),
        axis.title.y=element_blank()) +
  facet_grid(rows=vars(roc_df$feature_name),scales = "free_y",space='free_y' ,switch = 'y')+
  scale_x_discrete(expand = c(0,0)) -> pp
  #scale_y_discrete())
library(grid)

pdf(file='genomic_heatmap.pdf',width = 8.5,height = 11)
#grid::grid.newpage()
#vp <- viewport(width = unit(7.5, "inches"), height = unit(10, "inches"))
#pushViewport(vp)
grid::grid.draw(pp)
y1<-0.07
dy<-0.22
ddy<-0.006
dx <- 0.01
xx <- 0.15
ll <-1
grid.lines(x=c(xx,xx),y=c(y1,y1+dy),gp=grid::gpar(lwd=ll))
grid.lines(x=c(xx,xx+dx),y=c(y1,y1),gp=grid::gpar(lwd=ll))
grid.lines(x=c(xx,xx+dx),y=c(y1+dy,y1+dy),gp=grid::gpar(lwd=ll))
for (val in 1:3) {
  y1<-y1+ddy+dy
  grid.lines(x=c(xx,xx),y=c(y1,y1+dy),gp=grid::gpar(lwd=ll))
  grid.lines(x=c(xx,xx+dx),y=c(y1,y1),gp=grid::gpar(lwd=ll))
  grid.lines(x=c(xx,xx+dx),y=c(y1+dy,y1+dy),gp=grid::gpar(lwd=ll))
}
dev.off()

