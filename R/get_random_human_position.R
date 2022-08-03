#' Title
#'
#' @param nn number of randon site needed
#'
#' @return a tibble with the sites
#'
get_random_human_positions <- function(nn=10) {

  genome_sequence <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  genome_sequence@user_seqnames <- genome_sequence@user_seqnames[genome_sequence@user_seqnames %in% paste0("chr", c(1:22, "X", "Y", "M"))]
  genome_sequence@seqinfo <- genome_sequence@seqinfo[paste0("chr", c(1:22, "X", "Y", "M"))]


  set.seed(as.numeric(Sys.time()))

  chr<- sample(names(seqlengths(genome_sequence)),nn,replace = TRUE)
  strand <- sample(c('-','+'),nn,replace = TRUE)
  seqlengths(genome_sequence)[chr]
  position <- sapply(seqlengths(genome_sequence)[chr], function(xx){
      sample.int(xx,1)
  })

  tibble::tibble(
    seqnames=chr,
    start=position,
    end=position,
    strand=strand
  )

  # tttt <- ttt %>%
  #   dplyr::rename(start=.data$position) %>%
  #   dplyr::mutate(end=.data$start) %>%
  #   dplyr::mutate(siteID=NULL) %>%
  #   dplyr::relocate(.data$strand,.after = tidyselect::last_col()) %>%
  #   dplyr::mutate(patient=rep(group_vector,3)) %>%
  #   dplyr::mutate(type='match') %>%
  #   dplyr::rename(seqnames=.data$chr)
}
