% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/support.R
\name{point_ang_r}
\alias{point_ang_r}
\title{Find coordinates of point some path away from origin}
\usage{
point_ang_r(p, r, theta)
}
\arguments{
\item{p}{numeric, the origin. A vector or two values, x and y of a point in cartesian space.}

\item{r}{numeric, the path length or radius.}

\item{theta}{numeric, the angle of the path.}
}
\value{
A vector of length 2 represent x and y position
}
\description{
Internal function used to find x,y positions of paths given and origin, radius, and angle.
}
