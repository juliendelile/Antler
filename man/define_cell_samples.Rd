% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{define_cell_samples}
\alias{define_cell_samples}
\title{Define the cells sample names and map the cell colors to some of their metadata.}
\arguments{
\item{sample_feature}{a character vector of cell features. Must be included in the column names of the cell metadata matrix.}

\item{color_mapping}{a character vector having the same length as 'sample_features' and whose values are all different and either 'Hue', 'Saturation' or 'Value'. If NULL (default) the colors are chosen randomly.}

\item{seed}{an integer setting the seed for random color generation.}
}
\description{
Define the cells sample names and map the cell colors to some of their metadata.
}
