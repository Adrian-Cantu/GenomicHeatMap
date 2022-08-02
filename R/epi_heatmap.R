
#' Title of the song
#'
#' @return naive declaration
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom tibble data_frame
#'
#' @examples
#' x <- epi_heatmap()
epi_heatmap <- function(){

  cc <- c('DnaseUwJurkat','H3K79me2','PolII')

  intSites_coor <- ViiV_IntSites %>%
    BiocGenerics::as.data.frame() %>%
    dplyr::filter(to_heatmap) %>%
    dplyr::select(c(seqnames,start,end,strand,patient,GTSP,Drug,concentration_nM,concentration_txt,Replicate)) %>%
    dplyr::mutate(type='insertion') %>%
    dplyr::group_by(seqnames) %>%
    dplyr::slice_sample(n=30,replace = TRUE)

  #set how the groups are defined by changeing patient
  group_vector<-paste0(intSites_coor$Drug,'_',intSites_coor$concentration_txt)
  intSites_coor$patient <- group_vector


  df_to_randomize <- tibble::tibble(.rows = nrow(intSites_coor)) %>%
    dplyr::mutate(siteID=dplyr::row_number()) %>%
    dplyr::mutate(gender='m')


  genome_sequence <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  genome_sequence@user_seqnames <- genome_sequence@user_seqnames[genome_sequence@user_seqnames %in% paste0("chr", c(1:22, "X", "Y", "M"))]
  genome_sequence@seqinfo <- genome_sequence@seqinfo[paste0("chr", c(1:22, "X", "Y", "M"))]


  set.seed(as.numeric(Sys.time()))
  ## really need to remove this dependency
  ttt<-intSiteRetriever::get_N_MRCs(df_to_randomize,genome_sequence)

  tttt <- ttt %>%
    dplyr::rename(start=position) %>%
    dplyr::mutate(end=start) %>%
    dplyr::mutate(siteID=NULL) %>%
    dplyr::relocate(strand,.after = last_col()) %>%
    dplyr::mutate(patient=rep(group_vector,3)) %>%
    dplyr::mutate(type='match') %>%
    dplyr::rename(seqnames=chr)

  to_get_features <- rbind(intSites_coor %>% dplyr::select(colnames(tttt)) ,tttt)


  #epi_files <- paste0(cc,'.rds')
  #names(epi_files) <- cc

  return(to_get_features)

}
