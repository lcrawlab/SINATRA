% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gp_inference.R
\name{find_lasso_variables}
\alias{find_lasso_variables}
\title{Conducts the Lasso Feature Selection}
\usage{
find_lasso_variables(data, radius = 0, weights = FALSE)
}
\arguments{
\item{data}{(nxm matrix) A nxm matrix containing the covariates and responses.}

\item{radius}{(positive integer ). An integer parameter specifying if, and how many neighboring features should be considered for feature selection.
This is done to capture critical points that may be 'close' to the selected feature.}
}
\value{
The output is a vector of indices/features to be selected.
}
\description{
\code{find_lasso_variables} returns a vector of the selected features/indices of the curve

The function uses Lasso to conduct feature selection. A cross validated lasso model is fit on the binary response data, and the feature importances are extracted
As lasso sets many of the features to 0, the selected features are those that are non-zero.
These features are returned along with a window of neighboring features if the user wishes (typically recommended).
}
