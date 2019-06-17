\name{v.chaingraph}
\alias{v.chaingraph}
\title{
Node-wise selection for (Two-layered) Gaussian graphical model with no symmetric constraint on precision matrix
}
\description{
This runs for a single node
}
\usage{
v.chaingraph(v,chlist,palist,Y,eta.prob=0.1,gamma.prob=0.1,lambda,delta,burnin.S,inf.S)
}
\arguments{
  \item{v}{
an integer for a target indice for node-wise regression
}
  \item{chlist}{
list for node numbers for layers
}
  \item{palist}{
list for node number for parent layers
}
  \item{Y}{
nxp matrix that includes genomic measurements for p genes/proteins and n number of samples
}
  \item{eta.prob}{
P(eta=1)
}
  \item{gamma.prob}{
P(gamma=1)
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
\item{Gamma}{|v.ch| x |v.pa| x inf.S dimensional array for between layer structure}
\item{eta}{|v.ch| x |v.ch| x inf.S dimensional array for within layer structure}
\item{A}{|v.ch| x |v.ch| x inf.S dimensional array for regression coefficients corresponding to edges within layer}
\item{B}{|v.ch| x |v.pa| x inf.S dimensional array for regression coefficients corresponding to edges between layers}
\item{kappa}{inf.S x |v.ch| matrix for the inverse variance}
}





