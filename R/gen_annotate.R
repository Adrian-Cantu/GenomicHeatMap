#' add genomic annotations
#'
#' @param ff a genomic ranges object
#'
#' @param features genomic features to calculate, any mix of c('GC_content')
#'
#' @param windows size of the windows to calculate
#'
#' @return the same genomic ranges object with genomic features added as extra columns
#' @export
#'
gen_annotate_df <- function(ff,features=c('GC_content','CpG','DNaseI_count'),windows=c(1000,10000,1000000)){
  # TODO need to add code to select genome and seqinfo
  genome_sequence <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  #genome_sequence@user_seqnames <- genome_sequence@user_seqnames %>% keepStandardChromosomes(pruning.mode="coarse")
  #genome_sequence@seqinfo <- genome_sequence@seqinfo %>% GenomeInfoDb::keepStandardChromosomes(pruning.mode="coarse")
  ff_df<- ff %>% as.data.frame()
  if (class(ff)[1]=='GRanges') {
    ff_gr <- ff %>% GenomeInfoDb::keepStandardChromosomes(pruning.mode="coarse")
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
    if(xx=='GC_content'){
      for(nn in names(windows)) {
        gg_esp <- (ff_gr+windows[nn]/2) %>% GenomicRanges::trim()
        GenomeInfoDb::isCircular(gg_esp) <- rep(FALSE,length(GenomeInfoDb::isCircular(gg_esp)))
        gg_esp <- gg_esp %>% GenomicRanges::trim()
        tt_genome_seq <- genome_sequence
        # TODO this is a horrible hack
        si <- GenomeInfoDb::seqinfo(tt_genome_seq)
        GenomeInfoDb::isCircular(si)["chrM"] <- FALSE
        tt_genome_seq@seqinfo <- si
        BSgenome::getSeq(tt_genome_seq,gg_esp) %>%
          Biostrings::alphabetFrequency() %>%
          tibble::as_tibble() %>% #TODO return gc=0 for strings with no ACGTs
          dplyr::mutate(gc=(.data$C+.data$G)/(.data$A+.data$T+.data$C+.data$G)) %>%
          dplyr::pull(gc) ->
          gen_df[,paste('GC_content',nn,sep = '.')]
      }
    }
    if(xx=='CpG'){
      bsession <- rtracklayer::browserSession()
      rtracklayer::genome(bsession) <- "hg38"
      ttt <- rtracklayer::ucscTableQuery(bsession, track = "CpG Islands", table = "cpgIslandExt") %>%
        rtracklayer::getTable() %>%
        GenomicRanges::makeGRangesFromDataFrame() %>%
        GenomeInfoDb::keepStandardChromosomes(pruning.mode="coarse")
      for(nn in names(windows)) {
        GenomicRanges::countOverlaps(ff_gr,ttt, maxgap = windows[nn]/2, ignore.strand = TRUE) ->
          gen_df[,paste('CpG',nn,sep = '.')]
      }
    }
    if(xx=='DNaseI_count'){
      bsession <- rtracklayer::browserSession()
      rtracklayer::genome(bsession) <- "hg38"
      ttt <- rtracklayer::ucscTableQuery(bsession, track = "DNase Clusters", table = "wgEncodeRegDnaseClustered") %>%
        rtracklayer::getTable() %>%
        GenomicRanges::makeGRangesFromDataFrame() %>%
        GenomeInfoDb::keepStandardChromosomes(pruning.mode="coarse")
      for(nn in names(windows)) {
        GenomicRanges::countOverlaps(ff_gr,ttt, maxgap = windows[nn]/2, ignore.strand = TRUE) ->
          gen_df[,paste('DNaseI_count',nn,sep = '.')]
      }
    }
  }
  return(cbind(ff_df,gen_df) %>% GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE))
}
