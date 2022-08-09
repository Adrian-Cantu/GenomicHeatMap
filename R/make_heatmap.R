#' Title adds anotation columns to dataframe
#'
#' @param feature_df datafram with real and random sites
#'
#' @return tum
#' @export
#'
#' @importFrom GenomeInfoDb keepSeqlevels
#'
epi_annotate_df <- function(feature_df){

  epifiles_tmp <- list.files(system.file("exdata", package = "GenomicHeatMap"))
  epifiles <- epifiles_tmp[grepl('\\.rds$',epifiles_tmp)]
  names(epifiles) <- epifiles %>% stringr::str_remove(".rds")

  kk <- purrr::imap(epifiles, function(x,name){
    xx <- readRDS(system.file("exdata",x,package = "GenomicHeatMap"))
    feature_df <<- hiAnnotator::getFeatureCounts(feature_df, xx, name)
    rm(xx)
    gc()
    return(1)
  })
  return(feature_df)
}


#' Title sort feature
#'
#' @param ... a vector
#'
#' @return a sorted vector
#' @importFrom rlang .data
#' @importFrom magrittr %>%
sort_features <- function(...){
  translate_window <- stats::setNames(c(1000,1000000), c("Kb", "Mb"))
  xx <- tibble::tibble(.rows = length(unique(as.vector(...)))) %>%
    dplyr::mutate(site=unique(as.vector(...))) %>%
    tidyr::separate(.data$site,into=c('fname','window'),sep = '\\.',remove = FALSE) %>%
    tidyr::separate(.data$window,into=c('size','todel'),sep = '\\D+$',remove = FALSE) %>%
    tidyr::separate(.data$window,into=c('todel','size_word'),sep = '^\\d+',remove = FALSE) %>%
    dplyr::mutate(size_mult=translate_window[.data$size_word]) %>%
    dplyr::mutate(size_sort=as.numeric(.data$size)*.data$size_mult) %>%
    dplyr::mutate(todel=NULL) %>%
    dplyr::arrange(.data$fname,rev(.data$size_sort)) %>%
    dplyr::pull(var = .data$site)
  toret <- factor(as.vector(...),xx)
  return(toret)
}



#' Title makes roc object from dataframe
#'
#' @param annotate_df a dataframe
#'
#' @return ROC object
#' @export
#'
#'
#' @importFrom rlang .data
#' @importFrom magrittr %>%

make_roc <- function(annotate_df){
  #utils::globalVariables("where")
  to_roc_df <-as.data.frame(annotate_df) %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~tidyr::replace_na(.x, 0)))
  #%>% replace_na(0) #replace(is.na(.data), 0)

  roc.res <- hotROCs::ROC.ORC(
    response = to_roc_df$type,
    variables = to_roc_df %>% dplyr::select(-c(.data$seqnames,.data$start,.data$end,
                                        .data$width,.data$strand,.data$heat_group,.data$type)),
    origin=to_roc_df$heat_group)

  return(roc.res)
}


#' makes heatmap figures
#'
#' @param roc.res roc object
#' @param title title to display on the figure
#'
#' @return a ggplot
#' @export
#'
#' @importFrom rlang .data
make_heatmap <- function(roc.res,title='heatmap'){

  roc_df <- roc.res$ROC %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "feature") %>%
    tidyr::pivot_longer(!.data$feature,values_to='val', names_to='sample') %>%
    dplyr::mutate(feature=sort_features(.data$feature)) %>%
    tidyr::separate(.data$feature,into = c('feature_name','feature_concentration'),sep = '\\.',remove = FALSE)


  roc_df %>%
    ggplot2::ggplot( ggplot2::aes(y=.data$feature,x=.data$sample, fill= .data$val)) +
    ggplot2::geom_tile() +
    ggplot2::theme_classic() +
    ggplot2::labs(fill="ROC area",title=title) +
    ggplot2::scale_fill_gradientn(colours=c('blue','grey90','red'),
                         na.value = "transparent",
                         breaks=c(0,0.5,1),
                         labels=c(0,0.5,1),
                         limits=c(0,1)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30,vjust = 1, hjust=1),
          axis.text.y.left = ggplot2::element_text(size=9),
          axis.title.x=ggplot2::element_blank(),
          panel.spacing.y = grid::unit(-0.15, "line"),
          strip.placement='outside',
          panel.border = ggplot2::element_blank(),
          panel.background= ggplot2::element_blank(),
          strip.background = ggplot2::element_blank(),
          strip.text.y.left = ggplot2::element_text(angle=0,size=10),
          plot.background = ggplot2::element_blank(),
          axis.ticks.length.y.left=grid::unit(0.1,'line'),
          axis.title.y=ggplot2::element_blank()) +
    ggplot2::facet_grid(rows=ggplot2::vars(roc_df$feature_name),scales = "free_y",switch = 'y')+
    ggplot2::scale_x_discrete(expand = c(0,0)) +
    ggplot2::scale_y_discrete(labels=roc_df$feature_concentration,breaks=roc_df$feature,expand = c(0,0))

}
