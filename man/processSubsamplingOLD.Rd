% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrapperProcessData.R
\name{processSubsamplingOLD}
\alias{processSubsamplingOLD}
\title{Sub-sample and split data into training and test set}
\usage{
processSubsamplingOLD(
  cluster_selection_out,
  subsample_num = 1000,
  train_test_ratio = 0.5,
  cluster_proportion = "proportional",
  verbose = TRUE,
  seed = 8,
  ...
)
}
\arguments{
\item{cluster_selection_out}{Output of `processInputFormat`.}

\item{subsample_num}{Number of cells after sub-sammpling.}

\item{train_test_ratio}{Training to test data ratio.}

\item{cluster_proportion}{\itemize{
  \item \code{proportional:} (default) Same proportion of cells in each cluster, for the training and data sets, compared to the original cluster proportion.
  \item \code{equal}: Equal proportion of cells in each cluster, for the training and data sets, compared to the original cluster proportion.
}}
}
\value{
List containing
\itemize{
  \item \code{training_matrix}:
  \item \code{training_clusters}:
  \item \code{training_clusters_num}:
  \item \code{test_matrix}:
  \item \code{test_clusters}:
  \item \code{test_clusters_num}:
}
}
\description{
Sub-sample and split data into training and test set
}
