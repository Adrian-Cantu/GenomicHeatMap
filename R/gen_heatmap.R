# library
library(tidyverse)
library(hiAnnotator)
library(purrr)
library(hotROCs)
library(intSiteRetriever)
#get_N_MRCs()
library(GenomicRanges)


#
genome_sequence <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
genome_sequence@user_seqnames <- genome_sequence@user_seqnames[genome_sequence@user_seqnames %in% paste0("chr", c(1:22, "X", "Y", "M"))]
genome_sequence@seqinfo <- genome_sequence@seqinfo[paste0("chr", c(1:22, "X", "Y", "M"))]
#

intSites <- readRDS("20220708_ViiV_insites_plus_sampleinfo.rds")

### need to chose the right grouping variable, maybe patient or GTSP or some other
intSites_coor <- intSites %>%
  as.data.frame() %>%
  filter(to_heatmap) %>%
  select(c(seqnames,start,end,strand,patient,GTSP,Drug,concentration_nM,concentration_txt,Replicate)) %>%
  mutate(type='insertion')


### start roc code

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
window_size_refSeq <- c("10k"=1e4, "100k"=1e5, "1M"=1e6)
window_size_CpG_counts <- c("1k"=1e3, "10k"=1e4)
window_size_CpG_density <- c("10k"=1e4, "100k"=1e5, "1M"=1e6)
window_size_GC <- c("100"=100, "1k"=1000, "10k"=1e4, "100k"=1e5, "1M"=1e6)
#window_size_GC <- c("100"=100, "1k"=1000)
window_size_DNaseI <- c("1k"=1e3, "10k"=1e4, "100k"=1e5, "1M"=1e6)
window_size_epi <- c("10k"=1e4)

#
refGenes <- readRDS(file="hg38.refSeq.rds")
refGenes <- refGenes[seqnames(refGenes) %in% paste0("chr", c(1:22, "X", "Y", "M"))]

CpG_data <- cpg <- getUCSCtable("cpgIslandExt", "CpG Islands", freeze = "hg38") %>%
  dplyr::filter(chrom %in% paste0("chr", c(1:22, "X", "Y", "M")))

#getUCSCtable("cpgIslandExt", "CpG Islands", freeze = "hg38")


CpG_islands <- GenomicRanges::GRanges(
  seqnames = CpG_data$chrom,
  ranges = IRanges::IRanges(
    start = CpG_data$chromStart, end = CpG_data$chromEnd
  ),
  strand = "*",
  seqinfo = GenomeInfoDb::seqinfo(genome_sequence)
)

mcols(CpG_islands) <- CpG_data
#
DNaseI_data <- getUCSCtable("wgEncodeRegDnaseClustered", "DNase Clusters", freeze = "hg38") %>%
  dplyr::filter(chrom %in% paste0("chr", c(1:22, "X", "Y", "M")))

DNaseI <- GenomicRanges::GRanges(
  seqnames = DNaseI_data$chrom,
  ranges = IRanges::IRanges(start = DNaseI_data$chromStart, end = DNaseI_data$chromEnd),
  strand = "*",
  seqinfo = GenomeInfoDb::seqinfo(genome_sequence))

mcols(DNaseI) <- DNaseI_data
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
####

####
#all_names <- unique(intSites$GTSP)[100:104]
all_names <- unique(to_get_features$patient)
#plan(sequential)
#plan(multisession, workers = .num_cores)
l_names <- length(all_names)

hg38_seqinfo <- genome_sequence@seqinfo
#hg38_seqlev <- seqlevels(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)


  c_sample <- to_get_features

# 142938751

  gen_final_data <- c_sample %>% as.data.frame() %>% sample_n(4000) %>%   #filter(start==142938751) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=TRUE,seqinfo = hg38_seqinfo) %>%
    hiAnnotator::getFeatureCounts(refGenes, "refSeq_counts", width = window_size_refSeq) %>%
    GCcontent::getGCpercentage("GC", window_size_GC, genome_sequence) %>%
    hiAnnotator::getFeatureCounts(CpG_islands, "CpG_counts", width = window_size_CpG_counts) %>%
    hiAnnotator::getFeatureCounts(CpG_islands, "CpG_density", width = window_size_CpG_density) %>%
    from_counts_to_density("CpG_density", window_size_CpG_density) %>%
    hiAnnotator::getFeatureCounts(DNaseI, "DNaseI_count", width = window_size_DNaseI) %>%
    as.data.frame()

  hg38_seqinfo %>% as.data.frame()


gen_final_data %>% filter(is.na(DNaseI_count.1k)) %>% pull(var=seqnames) %>%  unique()
# chr7 chr8 chr3 chr4 chr6 chr2 chr9 chr5 chrX


gen_final_data %>% filter(if_any(everything(), is.na)) -> kk
roc.res <- ROC.ORC(
  response = gen_final_data$type,
  variables = gen_final_data %>% select(-c(seqnames,start,end,width,strand,patient,type)),
  origin=gen_final_data$patient)
