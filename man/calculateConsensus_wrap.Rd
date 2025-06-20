% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrapperMethods.R
\name{calculateConsensus_wrap}
\alias{calculateConsensus_wrap}
\title{Wrapper function to calculate consensus whilst selecting certain methods}
\usage{
calculateConsensus_wrap(
  list_markers_temp,
  final_out,
  num_markers = num_markers,
  chosen_measure = "F1_macro",
  list_performance_valid = list(),
  metric_thres = 0,
  metric_topnum = 1,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{list_markers_temp}{List of string vector, with each element of list corresponding to a marker method}

\item{final_out}{List of matrix and cluster information produced by function `processSubsampling`}

\item{num_markers}{Number of markers}

\item{chosen_measure}{The performance measure used to choose the methods used in the consensus. Options are precision_weighted, precision_macro, recall_weighted, recall_macro, F1_macro, F1_weighted, precision_micro}

\item{list_performance_valid}{Performance of the chosen markers using the validation matrix. Each element in list corresponds to a method.}

\item{metric_thres}{Threshold above which a method will be considered in the consensus calculation, based on `chosen_measure`}

\item{metric_topnum}{Number of methods to consider in consensus (based on sorted `chosen_measure`).}
}
\value{
Consensus markers
}
\description{
Wrapper function to calculate consensus whilst selecting certain methods
}
