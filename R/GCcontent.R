#' calulate GC for several window sizes
#'
#' @param reference_genome_sequence actual reference genome (from BSgenome)
#' to get reference_genome_sequence use get_reference_genome() function from
#' random_site.R
#'
#' @note start up time to load human genome is about 10 second
getGCpercentage <- function(
    sites, column_prefix, window_size, reference_genome_sequence
){
  stopifnot(length(window_size) == length(names(window_size)))
  metadata <-  GenomicRanges::mcols(sites)

  rangesToCalc <- expand_trim_GRanges(sites, window_size)

  #seqs will take a lot of memory
  #could split at severe cpu time penelty
  seqs <- getSeq(reference_genome_sequence, rangesToCalc, as.character=F)

  letterFreqs <- letterFrequency(seqs, c("G", "C", "A", "T"))
  rm(seqs)

  GC <- letterFreqs[, c("G", "C")]
  ATGC <- letterFreqs[, c("A", "T", "G", "C")]

  gcContent <- rowSums(GC)/rowSums(ATGC)

  gcContent[!is.finite(gcContent)] <- NA #handled gracefully by pipeUtils

  gcContent <- DataFrame(matrix(gcContent, nrow=length(sites)))

  names(gcContent) <- paste(column_prefix, names(window_size), sep=".")

  return(gcContent)
  # mcols(sites) <- cbind(metadata, gcContent)
  #
  # sites
}

expand_trim_GRanges <- function(sites, window_size) {
  nsites <- length(sites)
  strand(sites) = "+" #unimportant for GC and speeds up later calculations

  sites.seqinfo.original <- seqinfo(sites)
  isCircular(seqinfo(sites)) <- rep(FALSE, length(seqinfo(sites)))

  sites <- rep(sites, length(window_size))
  sites <- trim(suppressWarnings(flank(sites,
                                       rep(window_size/2, each=nsites),
                                       both=T)))

  seqinfo(sites) = sites.seqinfo.original
  sites
}
