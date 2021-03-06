% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulated_data_generation.R
\name{generate_normal_complex_fixed}
\alias{generate_normal_complex_fixed}
\title{Generate Kernel interpolated meshes.}
\usage{
generate_normal_complex_fixed(
  grid_size = 25,
  noise_points,
  causal_points1,
  causal_points2,
  func = rbf_gauss,
  eta = 5
)
}
\arguments{
\item{grid_size}{(int) : The fineness of the grid for which the interpolation is on.}

\item{noise_points}{(matrix) : Matrix of points to be shared across both classees of shapes.}

\item{causal_points1}{(matrix) : Matrix of points unique to class 1.}

\item{causal_points2}{(matrix) : Matrix of points unique to class 2.}

\item{func}{(function)  : The function used for the kernel interpolation.}

\item{eta}{(float)  : The shape paramter for the kernel.}
}
\value{
complex (list) : List of data containing the mesh versions of the interpolated shapes, causal points, and the shared points.
}
\description{
Given a set of input data of shared points and causal points,
 the function generates two classes of interpolated 3d shapes with a specified kernel.
}
