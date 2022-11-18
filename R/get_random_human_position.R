#' Title
#'
#' @param nn number of randon site needed
#' @param r_seed seed to use for generation of random control insertion sites
#'
#' @return a tibble with the sites
#'
get_random_human_positions <- function(nn=10,r_seed=as.numeric(Sys.time())) {

  genome_sequence <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  genome_sequence@user_seqnames <- genome_sequence@user_seqnames[genome_sequence@user_seqnames %in% paste0("chr", c(1:22, "X", "Y", "M"))]
  genome_sequence@seqinfo <- genome_sequence@seqinfo[paste0("chr", c(1:22, "X", "Y", "M"))]


  set.seed(r_seed)

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
}

#' Get a random position in the mac5 genome
#'
#' @param nn number of randon site needed
#' @param r_seed seed to use for generation of random control insertion sites
#'
#' @return a tibble with the sites
#'
get_random_macFas5_positions <- function(nn=10,r_seed=as.numeric(Sys.time())) {

  genome_sequence <- BSgenome.Mfascicularis.NCBI.5.0::BSgenome.Mfascicularis.NCBI.5.0
  genome_sequence@user_seqnames <- genome_sequence@user_seqnames %>% stringr::str_replace('^MFA','chr') %>% stringr::str_replace('^MT','chrM')
  all_seqlenghts <- seqlengths(genome_sequence)[names(seqlengths(genome_sequence)) %in% paste0("chr", c(1:20, "X", "Y", "M"))]
  genome_sequence@user_seqnames <- genome_sequence@user_seqnames[genome_sequence@user_seqnames %in% paste0("chr", c(1:20, "X", "Y", "M"))]
  genome_sequence@seqinfo <- genome_sequence@seqinfo[paste0("chr", c(1:20, "X", "Y", "M"))]


  set.seed(r_seed)

  chr<- sample(names(seqlengths(genome_sequence)),nn,replace = TRUE)
  strand <- sample(c('-','+'),nn,replace = TRUE)
  #seqlengths(genome_sequence)[chr]
  position <- sapply(all_seqlenghts[chr], function(xx){
    sample.int(xx,1)
  })

  tibble::tibble(
    seqnames=chr,
    start=position,
    end=position,
    strand=strand
  )
}



#' Get a random position in an specified genome genome
#'
#' @param nn number of randon site needed
#' @param r_seed seed to use for generation of random control insertion sites
#' @param genome genome to get random sites from. one of c('hg38','macFas5')
#'
#' @return a tibble with the sites
#'
get_random_positions <- function(nn=10,r_seed=as.numeric(Sys.time()),genome='hg38') {
  if (genome=='hg38'){
    get_random_human_positions(nn=nn,r_seed = r_seed)
  } else if (genome=='macFas5') {
    get_random_macFas5_positions(nn=nn,r_seed = r_seed)
  } else{
    get_random_macFas5_positions(nn=nn,r_seed = r_seed)
  }

}
