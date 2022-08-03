
#' Title of the song
#'
#' @return naive declaration
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom tibble data_frame
#' @importFrom rlang .data
#' @examples
#' x <- epi_heatmap()
epi_heatmap <- function(){

  cc <- c('DnaseUwJurkat','H3K79me2','PolII')

  intSites_coor <- ViiV_IntSites %>%
    BiocGenerics::as.data.frame() %>%
    dplyr::filter(.data$to_heatmap) %>%
    dplyr::select(c(.data$seqnames,.data$start,.data$end,.data$strand,.data$GTSP,.data$Drug,
                    .data$concentration_nM,.data$concentration_txt,.data$Replicate)) %>%
    dplyr::mutate(type='insertion') %>%
    dplyr::mutate(patient=paste0(.data$Drug,'_',.data$concentration_txt)) %>%
    dplyr::group_by(.data$seqnames) %>%
    dplyr::slice_sample(n=30,replace = TRUE)

  #set how the groups are defined by changeing patient
  group_vector <- intSites_coor$patient

#
#   df_to_randomize <- tibble::tibble(.rows = nrow(intSites_coor)) %>%
#     dplyr::mutate(siteID=dplyr::row_number()) %>%
#     dplyr::mutate(gender='m')



  tttt <- get_random_human_positions( nn = intSites_coor %>% nrow() *3) %>%
    dplyr::mutate(patient=rep(group_vector,3)) %>%
    dplyr::mutate(type='match')

  to_get_features <- rbind(intSites_coor %>% dplyr::select(colnames(tttt)) ,tttt)


  #epi_files <- paste0(cc,'.rds')
  #names(epi_files) <- cc

  return(to_get_features)

}
