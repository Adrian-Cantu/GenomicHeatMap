
#' Title of the song
#'
#' @return naive declaration
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom rlang .data
epi_heatmap_test <- function(){

  cc <- c('DnaseUwJurkat','H3K79me2','PolII')

  intSites_coor <- ViiV_IntSites %>%
    BiocGenerics::as.data.frame() %>%
    dplyr::filter(.data$to_heatmap) %>%
    dplyr::select(c(.data$seqnames,.data$start,.data$end,.data$strand,.data$GTSP,.data$Drug,
                    .data$concentration_nM,.data$concentration_txt,.data$Replicate)) %>%
    dplyr::mutate(type='insertion') %>%
    dplyr::mutate(patient=paste0(.data$Drug,'_',.data$concentration_txt)) %>%
    dplyr::group_by(.data$seqnames) %>%
    dplyr::slice_sample(n=30,replace = TRUE) %>%
    dplyr::rename(heat_group=.data$patient)

  group_vector <- intSites_coor$heat_group
#
#   df_to_randomize <- tibble::tibble(.rows = nrow(intSites_coor)) %>%
#     dplyr::mutate(siteID=dplyr::row_number()) %>%
#     dplyr::mutate(gender='m')



  tttt <- get_random_human_positions( nn = intSites_coor %>% nrow() *3) %>%
    dplyr::mutate(heat_group=rep(group_vector,3)) %>%
    dplyr::mutate(type='match')

  to_get_features <- rbind(intSites_coor %>% dplyr::select(colnames(tttt)) ,tttt)


  #epi_files <- paste0(cc,'.rds')
  #names(epi_files) <- cc

  return(to_get_features)

}

#' Format
#'
#' @param intSites insertion site dataframe
#' @param group_variable variable that will be used to group the data, each value corresponds to a
#' column in the final heatmap
#' @param r_seed seed to use for generation of random control insertion sites
#' @param genome genome to get random sites from. one of c('hg38','macFas5')
#'
#' @return naive declaration
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom rlang .data
intsite_to_heatmap_df <- function(intSites,group_variable='patient',r_seed=as.numeric(Sys.time()),genome='hg38'){

  intSites_coor <- intSites %>%
    BiocGenerics::as.data.frame() %>%
    dplyr::select(c(.data$seqnames,.data$start,.data$end,.data$strand,.data[[group_variable]])) %>%
    dplyr::mutate(type='insertion') %>%
    dplyr::rename(heat_group=.data[[group_variable]]) #%>%
    #dplyr::filter(seqnames!='chrM')

  group_vector <- intSites_coor$heat_group

  tttt <- get_random_positions( nn = intSites_coor %>% nrow() *3,r_seed=r_seed,genome = genome) %>%
    dplyr::mutate(heat_group=rep(group_vector,3)) %>%
    dplyr::mutate(type='match')

  to_get_features <- rbind(intSites_coor %>% dplyr::select(colnames(tttt)) ,tttt)

  return(to_get_features)
}

