\name{ch.chaingraph.str}
\alias{ch.chaingraph.str}
\title{
Structured estimation for (Two-layered) Gaussian graphical model
}
\description{
Structured MCMC jointly for a layer
}
\usage{
ch.chaingraph.str(v.ch,v.pa,Y,G,lambda,delta,burnin.S,inf.S)
}
\arguments{
  \item{v.ch}{
a vector for indices of the target layer
}
  \item{v.pa}{
a vector for indices of the parent layer
}
  \item{Y}{
nxp matrix that includes genomic measurements for p genes/proteins and n number of samples
}
  \item{G}{
Adjacency matrix for the mlGGM across p variables
}
  \item{lambda}{
a integer value for hyperparameter for variance and undirected edges
}
  \item{delta}{
a integer value for hyperparameter for variance and undirected edges
}
  \item{burnin.S}{
length of MCMC for burn-in
}
  \item{inf.S}{
length of MCMC for inference
}
}

\value{
a list
\item{A}{|v.ch| x |v.ch| x inf.S dimensional array for regression coefficients corresponding to edges within layer}
\item{B}{|v.ch| x |v.pa| x inf.S dimensional array for regression coefficients corresponding to edges between layers}
\item{kappa}{inf.S x |v.ch| matrix for the inverse variance}
}





