#' Title adds anotation columns to dataframe
#'
#' @param feature_df dataframe with real and random sites
#' @param epi_directory directory with .rds of the tracks
#'
#' @return tum
#' @export
#'
#' @importFrom GenomeInfoDb keepSeqlevels
#'
epi_annotate_df <- function(feature_df,epi_directory = system.file("exdata", package = "GenomicHeatMap")){

  epifiles_tmp <- list.files(epi_directory)
  epifiles <- epifiles_tmp[grepl('\\.rds$',epifiles_tmp)]
  names(epifiles) <- epifiles %>% stringr::str_remove(".rds")

  kk <- purrr::imap(epifiles, function(x,name){
    xx <- readRDS(file = file.path(epi_directory,x))
    feature_df <<- hiAnnotator::getFeatureCounts(feature_df, xx, name)
    rm(xx)
    gc()
    return(1)
  })
  return(feature_df)
}


#' Title adds anotation columns to dataframe
#'
#' @param ff_df datafram with real and random sites
#' @param epi_directory directory with .rds of the tracks
#'
#' @return tum
#' @export
#'
#' @importFrom GenomeInfoDb keepSeqlevels
#'
epi_annotate_df2 <- function(ff_df,epi_directory = system.file("exdata", package = "GenomicHeatMap")){
  epifiles_tmp <- list.files(epi_directory)
  epifiles <- epifiles_tmp[grepl('\\.rds$',epifiles_tmp)]
  names(epifiles) <- epifiles %>% stringr::str_remove(".rds")
  epinames <- names(epifiles)
  ll_names <- outer(epinames, c('1Kb','10Kb','1Mb'), paste, sep=".") %>%
    t() %>%
    as.vector()
  epi_df <- matrix(0,nrow = nrow(ff_df),ncol = length(ll_names),dimnames = list(NULL,ll_names)) %>%
    tibble::as_tibble()
  ff_gr <- ff_df %>% GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

  for (ll in 1:length(epifiles)) {
    x <- epifiles[ll]
    name<-names(epifiles)[ll]
    xx <- readRDS(file=file.path(epi_directory,x))
    nums   <- c(1000,10000,1000000)
    nums_l <- c('1Kb','10Kb','1Mb')
    for (lll in 1:length(nums)) {
      GenomicRanges::countOverlaps(ff_gr,xx, maxgap = nums[lll]/2, ignore.strand = TRUE) ->
        epi_df[,paste(name,nums_l[lll],sep = '.')] #

    }
  }
  return(cbind(ff_df,epi_df) %>% GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE))
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
#'
#' @return a ggplot
#' @export
#'
#' @importFrom rlang .data
make_heatmap <- function(roc.res){

  roc_df <- roc.res$ROC %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "feature") %>%
    tidyr::pivot_longer(!.data$feature,values_to='val', names_to='sample') %>%
    dplyr::mutate(feature=sort_features(.data$feature)) %>%
    tidyr::separate(.data$feature,into = c('feature_name','feature_concentration'),sep = '\\.',remove = FALSE)

  # get stars for vehicle... code copied from Viiv project, needs much work -----------

  mk.stars <- function(pvmat) {
    x <- array("", dim(pvmat))
    x[] <- as.character(cut(pvmat, c(0, 0.001, 0.01, 0.05,
                                     1), c("***", "**", "*", ""), include.lowest = TRUE))
    x
  }

  roc.rows <- nrow(roc.res$ROC)
  roc.cols <- ncol(roc.res$ROC)

  opvals <- roc.res$pvalues$op
  opstars <- mk.stars(opvals)

  matchCol <- function(x, indx) x[2:1, ][x == indx]
  isCol <- function(x, indx) colSums(x == indx) == 1

  omasks <- lapply(1:roc.cols, function(x) {
    res <- array("", dim(roc.res$ROC))
    res[, x] <- "--"
    wc <- attr(opvals, "whichCol")
    res[, matchCol(wc, x)] <- opstars[, isCol(wc, x)]
    res
  })



  vehicle_stars <- omasks[[7]] %>%
    `colnames<-`(colnames(roc.res$ROC)) %>%
    `rownames<-`(rownames(roc.res$ROC)) %>%
    as.data.frame() %>%
    rownames_to_column(var = "feature") %>%
    pivot_longer(!feature,values_to='v_star', names_to='sample')

  # # to add stars as text to dataframe ---------------
  # ff_roc_df <- left_join(roc_df %>%  mutate(jj=paste0(feature,'_',sample)),
  #                        vehicle_stars %>% mutate(jj=paste0(feature,'_',sample)) %>% select(c(v_star,jj)),
  #                        by='jj') %>%
  #   mutate(jj=NULL)

  #-----------------------------------
  # code to get vs random stars

  #-----------------------------------
  stat_cuts <- c(0, 0.001, 0.01, 0.05, 1)
  roc_pval <- roc.res$pvalues$np %>%
    as.data.frame() %>%
    rownames_to_column(var = "feature") %>%
    pivot_longer(!feature,values_to='pval', names_to='sample') %>%
    mutate(feature=sort_features(feature)) %>%
    mutate(v_star=cut(pval, stat_cuts, labels = c("***", " **", " * ", "   "),include.lowest = TRUE))

  ff_roc_df <- left_join(roc_df %>% mutate(jj=paste0(feature,'_',sample)),
                          roc_pval%>% mutate(jj=paste0(feature,'_',sample)) %>% select(c(pval,v_star,jj)),
                          by='jj') %>%
    mutate(jj=NULL)



    return(ff_roc_df)
}

#' Title
#'
#' @param roc_df A dataframe generated from make_heatmap
#' @param title A title for the plot
#' @param star Whether to add pvalue significance as stars.
#'
#' @return a ggplot object
#' @export
#'
plot_heatmap <- function(roc_df,title='heatmap',star=FALSE) {

  plot_fig <- roc_df %>%
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

  if(star) {
    plot_fig <- plot_fig +
      ggplot2::geom_text(aes(label = v_star), color = "black", size = 3, nudge_y = -0.15)
  }

  return(plot_fig)

}
