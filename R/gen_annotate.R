#' add genomic annotations
#'
#' @param ff_df a genomic ranges
#'
#' @return the same genomic ranges object with genomic features added as extra columns
#' @export
#'
gen_annotate_df <- function(ff,features=c('GC_content'),windows=c(1000,10000,1000000)){
  # epifiles_tmp <- list.files(epi_directory)
  # epifiles <- epifiles_tmp[grepl('\\.rds$',epifiles_tmp)]
  # names(epifiles) <- epifiles %>% stringr::str_remove(".rds")
  # epinames <- names(epifiles)
  # TODO need to add code to select genome and seqinfo
  genome_sequence <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  #genome_sequence@user_seqnames <- genome_sequence@user_seqnames %>% keepStandardChromosomes(pruning.mode="coarse")
  genome_sequence@seqinfo <- genome_sequence@seqinfo %>% GenomeInfoDb::keepStandardChromosomes(pruning.mode="coarse")
  ff_df<- ff %>% as.data.frame()
  if (class(ff)[1]=='GRanges') {
    ff_gr <- ff %>% keepStandardChromosomes(pruning.mode="coarse")
    GenomeInfoDb::seqinfo(ff_gr) <- genome_sequence@seqinfo
  } else {
    ff_gr<- ff_df %>% GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE,seqinfo = genome_sequence@seqinfo)
  }
  names(windows) <-
    gdata::humanReadable(windows,standard = 'SI',sep='',justify = 'none',digits = 0) %>%
    stringr::str_replace('B$','')
  ll_names <- outer(features, names(windows), paste, sep=".") %>%
    t() %>%
    as.vector()
  gen_df <- matrix(0,nrow = nrow(ff_df),ncol = length(ll_names),dimnames = list(NULL,ll_names)) %>%
    tibble::as_tibble()
  # TODO this might break if length is < 1
  for (ll in 1:length(features)) {
    xx <- features[ll]
    #name<-names(epifiles)[ll]
    #xx <- readRDS(file=file.path(epi_directory,x))
    # for (lll in 1:length(nums)) {
    #   GenomicRanges::countOverlaps(ff_gr,xx, maxgap = nums[lll]/2, ignore.strand = TRUE) ->
    #     epi_df[,paste(name,nums_l[lll],sep = '.')] #
    #
    # }
    if(xx=='GC_content'){
      #window_size_GC <- c("1k"=1000, "10k"=1e4, "1M"=1e6)
      #gc_tt <- getGCpercentage("GC", window_size_GC, genome_sequence)
      #gen_df[,colnames(gc_tt)] <- gc_tt
      for(nn in names(windows)) {
        gg_esp <- (ff_gr+windows[nn]/2) %>% GenomicRanges::trim()
        #GenomeInfoDb::isCircular(gg_esp) <- rep(FALSE,length(GenomeInfoDb::isCircular(gg_esp)))
        # gg_esp <- gg_esp %>%
        #   GenomicRanges::trim()
        BSgenome::getSeq(genome_sequence,gg_esp) %>%
          Biostrings::alphabetFrequency() %>%
          tibble::as_tibble() %>% #TODO return gc=0 for strings with no ACGTs
          dplyr::mutate(gc=(.data$C+.data$G)/(.data$A+.data$T+.data$C+.data$G)) %>%
          dplyr::pull(gc) ->
          gen_df[,paste('GC_content',nn,sep = '.')]
      }
    }
  }
  return(cbind(ff_df,gen_df) %>% GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE))
}
