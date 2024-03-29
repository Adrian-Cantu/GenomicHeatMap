% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gen_annotate.R
\name{gen_annotate_df}
\alias{gen_annotate_df}
\title{add genomic annotations}
\usage{
gen_annotate_df(
  ff,
  features = c("GC_content", "CpG", "DNaseI_count"),
  windows = c(1000, 10000, 1e+06)
)
}
\arguments{
\item{ff}{a genomic ranges object}

\item{features}{genomic features to calculate, any mix of c('GC_content')}

\item{windows}{size of the windows to calculate}
}
\value{
the same genomic ranges object with genomic features added as extra columns
}
\description{
add genomic annotations
}
