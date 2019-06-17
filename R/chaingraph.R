giveA <- function(Afull,which) {
  # For undirected part
  if (length(which)>0) {
     matinv(Afull,which=which,negate=F)[-which,-which,drop=F]
	}else {NULL}
}

giveCov<- function(X,DD,kappa) {
# For undirected part
	p = ncol(X)
	n =nrow(X)
	if (p>0) {
		(diag(n) + tcrossprodCpp(prodCpp(X,diag(DD,nrow=p)),X))/kappa
	}else {diag(n)/kappa}
}

giveCov.dir <- function(X,CC,kappa) {
# For directed part
# Input:
#   - X:
#   - CCinv: inverse covariance for prior
#   - kappa: a scalar for inverse variance of the response
	p = ncol(X)
	n =nrow(X)
	if (p>0) {
		diag(n)/kappa +  tcrossprodCpp(prodCpp(X,diag(CC,nrow=p)),X)
	}else {diag(n)/kappa}
}

updateUndirected <- function(v, dat,lambda,delta,Alpha,eta,kappa,pmat,no.de,bCb,no.tau) {
  # sampling alpha and kappa given eta
  # Input
  #   - v: current node
  #   - dat: nxp data matrix
  #   - Afull: chol2inv(chol((crossprod(dat,tdat)+diag(p)*(1/lambda))))
  #   - lambda, delta: hyperparameters for kappa
  #   - eta: pxp eta matrix
  #   - kappa: px1 kappa vector (inverse covariance for all nodes)
  #   - pmat: pxp matrix for hyper parameters for eta
  #   - no.de : number of directed edges connected to each node (for updating kappa) (vector)
  #   - no.de : t(b) %*% inv(C) %*% b for updating kappa for each node (vector)
  #   - no.tau:  number of nodes in the current layer

  n = nrow(dat)
  p = ncol(dat)
  ### Graph movement ###
  if (sum(eta)==0) {is.swap=0
  }else{is.swap = rbinom(1,1,1/2)}
  is.move = FALSE
  if (!is.swap == 1) { #### ADD-DELETE
    w.ad = sample((1:p)[-v],1)
    new.eta = eta
    new.eta[v,w.ad] =new.eta[w.ad,v]= 1-eta[v,w.ad]
    w.up.l =  c(v,w.ad)## nodes for likelihood update
    is.move = TRUE
  }else{ #### SWAP
    cand0 = setdiff(which(eta[v,]==0),v)
    cand1 = which(eta[v,]==1)
    if (length(cand0)>0&length(cand1)>0) {
      w.ad1 = cand0[sample.int(length(cand0))[1]]
      w.ad2 = cand1[sample.int(length(cand1))[1]]
      new.eta = eta
      new.eta[v,w.ad1] = new.eta[w.ad1,v] = 1
      new.eta[v,w.ad2] = new.eta[w.ad2,v] = 0
      w.up.l = c(v,w.ad1,w.ad2)
      is.move = TRUE
    }
  }

  ### Sampling parameters ###
  if (is.move) {

	Afull = chol2inv(chol((crossprodCpp(dat,dat)+diag(p)*(1/lambda))))

	new.addr = sapply(w.up.l,function(x) which(new.eta[x,]==1),simplify=F)
	addr = sapply(w.up.l,function(x)which(eta[x,]==1),simplify=F)

	new.addr.inv = sapply(w.up.l,function(x) which(new.eta[x,]==0),simplify=F)
	addr.inv = sapply(w.up.l,function(x)which(eta[x,]==0),simplify=F)

	A.new = lapply(new.addr.inv,function(x) giveA(Afull,which=x))
	A =  lapply(addr.inv,function(x) giveA(Afull,which=x))

    #sample eta
	  lhr1 = sum(sapply(1:length(w.up.l),function(w) dmvnrm_arma(x = matrix(dat[,w.up.l[w]],ncol=n),mean=rep(0,n)
				,sigma = giveCov(dat[,new.addr[[w]],drop=F],rep(lambda,length(new.addr[[w]])),kappa[w.up.l[w]]),log=TRUE)))
	  hyper1 = sum(sapply(w.up.l,function(w) sum(new.eta[w,-w]*log(pmat[w,-w]) + (1-new.eta[w,-w])*log(1-pmat[w,-w]))))
	  lhr2 = sum(sapply(1:length(w.up.l),function(w) dmvnrm_arma(x = matrix(dat[,w.up.l[w]],ncol=n),mean=rep(0,n),
				sigma =giveCov(dat[,addr[[w]],drop=F],rep(lambda,length(addr[[w]])),kappa[w.up.l[w]]) ,log=TRUE)))
	  hyper2 = sum(sapply(w.up.l,function(w) sum(eta[w,-w]*log(pmat[w,-w]) + (1-eta[w,-w])*log(1-pmat[w,-w]))))
	  lhr = lhr1+hyper1-lhr2-hyper2

    if (log(runif(1))<lhr) {eta <- new.eta; A<-A.new}


    for (i in 1:length(w.up.l)) {
      w = w.up.l[i]
      alpha = rep(0,p)
      resid = dat[,w]
      if (!is.null(A[[i]])){
        # sample alpha
		addr = which(eta[w,]==1)
		  if (length(addr)>0) {
			  alpha[addr] = rmvnrm_arma(1,prodCpp(A[[i]],crossprodCpp(dat[,addr,drop=F],dat[,w,drop=F])),A[[i]]/kappa[w])
		  }
		Alpha[w,] = alpha
        # sample kappa
        #resid = dat[,w] -prodCpp(dat,alpha) Changed
        resid = dat[,w] -prodCpp(dat,as.matrix(alpha))
      }
      kappa[w] = rgamma(1,shape=(n+delta+no.tau-1+no.de[w]+sum(alpha!=0))/2,rate=(lambda +sum(resid^2)+bCb[w]+(lambda +no.tau-1)*sum(alpha^2) )/2)
    }
  } #if (is.move)
  return(list(eta=eta,kappa=kappa,Alpha=Alpha,is.move=is.move))
}

updateDirected <- function(y,X,lambda,delta,CCinv,B,Gamma,kappa,qmat,v.no.ue,v.aa,no.tau) {
  # sampling beta and kappa given gamma
  # Input
  #   - y: nx1 response vector corresponding to node v
  #   - X: nx(p*pP) covariate vector including parents of v and a function of parents of C(v)
  #   - lambda, delta: hyperparameters for kappa
  #   - CCinv: pxpP hyper parameter matrix for B
  #   - B : pxpP current B (permuted for response y)
  #   - Gamma: pxpP current gamma (permuted for response y)
  #   - kappa: a scalar for kappa for v (inverse covariance)
  #   - qmat: hyper parameter matrix for Gamma
  #   - v.no.ue: number of undirected edges for v (scalar) (for sampling kappa)
  #   - v.aa: t(alpha) %*% alpha (scalar) (for sampling kappa)
  #   - no.tau:  number of nodes in the current layer
    n = length(y)
    pP = ncol(B)
    p = nrow(B)
    vecpP = prod(dim(B))
    vqmat = c(qmat)
    ### Graph movement ###
    if (sum(B[1,])==0) {is.swap=0
    }else{is.swap = rbinom(1,1,1/2)}
    is.move = FALSE
    if (!is.swap == 1) { #### ADD-DELETE
      w.ad = sample(1:pP,1)
      new.Gamma = Gamma
      new.Gamma[1,w.ad] = 1-new.Gamma[1,w.ad]
      is.move = TRUE
    }else{ #### SWAP
      cand0 = which(Gamma[1,]==0)
      cand1 = which(Gamma[1,]==1)
      if (length(cand0)>0&length(cand1)>0) {
        w.ad1 = cand0[sample.int(length(cand0))[1]]
        w.ad2 = cand1[sample.int(length(cand1))[1]]
        new.Gamma = Gamma
        new.Gamma[1,w.ad1]  = 1
        new.Gamma[1,w.ad2]  = 0
        is.move = TRUE
      }
    }

    ### Sampling parameters ###
    if (is.move){
      vB = c(t(B))
      vGamma = c(t(Gamma))
      new.vGamma = c(t(new.Gamma))
	  CC = 1/diag(CCinv)
	  new.addr = which(new.vGamma==1)
	  addr = which(vGamma==1)
	  new.addr.inv = which(new.vGamma==0)
	  addr.inv = which(vGamma==0)

		if (pP*p>n) {
		  XC =  X%*%diag(CC)
		  tmp1 =chol2inv(chol(diag(n) + kappa * tcrossprodCpp(XC,X)))
		  Afull = diag(CC)-kappa*crossprodCpp(XC,prodCpp(tmp1,XC))
		}else{Afull = chol2inv(chol(kappa*crossprodCpp(X,X)+CCinv))}

		A.new = giveA(Afull,which=new.addr.inv)
		A = giveA(Afull,which=addr.inv)

      # sample gamma
		lhr1 = dmvnrm_arma(x = matrix(y,ncol=n),mean=rep(0,n)
							   ,sigma = giveCov.dir(X=X[,new.addr,drop=F],CC=CC[new.addr],kappa=kappa),log=TRUE)
		hyper1 = sum(new.vGamma * log(vqmat) + (1-new.vGamma) * log(1-vqmat))
		lhr2 = dmvnrm_arma(x=matrix(y,ncol=n),mean=rep(0,n),sigma=giveCov.dir(X=X[,addr,drop=F],CC=CC[addr],kappa=kappa),log=TRUE)
		hyper2 = sum(vGamma * log(vqmat) + (1-vGamma) * log(1-vqmat))

		lhr = lhr1+hyper1-lhr2-hyper2

      if (log(runif(1))<lhr) {vGamma <- new.vGamma; A<-A.new}

      # sample B
      vB = rep(0,length(vGamma))
      resid = y
      if (!is.null(A)){
		ww = which(vGamma==1)
		  if (length(ww)>0) {
				vB[ww]= rmvnrm_arma(1,kappa * prodCpp(A,crossprodCpp(X[,ww,drop=F],as.matrix(y))),A)
		  }
		#resid = y-prodCpp(X,vB)

        resid = y-prodCpp(X,as.matrix(vB))
      }
      kappa = rgamma(1,shape=(n+delta+no.tau-1+v.no.ue+sum(vB!=0))/2,rate=(lambda+sum(resid^2)+sum(vB^2/diag(CCinv))+(lambda +no.tau-1)*v.aa )/2)
      B = matrix(vB,ncol=pP,byrow=T)
      Gamma = matrix(vGamma,ncol=pP,byrow=T)
    }
    return(list(Gamma=Gamma,kappa=kappa,B=B,is.move=is.move))
}

ch.chaingraph <- function(v.ch,v.pa,Y,eta.prob=0.1,gamma.prob=0.1,lambda,delta,burnin.S,inf.S) {
    # This function fits regression model for a chain component
    # Input
    # - v.ch : indices for the target chain component
    # - v.pa : indices for parents set
    # - Y : nxp data matrix
    stopifnot(length(v.ch)>2)
    S = burnin.S + inf.S
    dat.C = scale(Y[,v.ch,drop=F])
    n = nrow(dat.C)
    p = ncol(dat.C)
    pmat = matrix(0.2,ncol(dat.C),ncol(dat.C))
    diag(pmat) = 0

    if (is.null(v.pa)) {
        Gamma.list = NULL
        B.list=NULL
    }else{
        dat.P = scale(Y[,v.pa,drop=F])
        pP = ncol(dat.P)
        Gamma.list = B.list = array(0,dim=c(p,pP,inf.S))
        CC = matrix(1/lambda,ncol=pP,nrow=p)# hyper parameter for b for nonzero gamma
        qmat = matrix(0.1,ncol(dat.C),ncol(dat.P))
    }
    eta.list =A.list= array(0,dim=c(p,p,inf.S))
    kappa.list = matrix(0,nrow=inf.S,ncol=p)
    ###########################################
    ####### Initialize all parameters #########
    ###########################################
    # eta (indicators for alpha) #
    w.upper = which(upper.tri(diag(p)))
    eta = matrix(0,p,p)
    eta[w.upper] = rbinom(length(w.upper),size=1,prob=eta.prob)
    eta = eta + tCpp(eta)
    diag(eta) = 0
    # Gamma (indicators for b)#
    if (!is.null(v.pa)) {Gamma = matrix(rbinom(p*pP,size=1,prob=gamma.prob),nrow=p,ncol=pP)}
    # A (pxp alpha)$
    A = pmat*eta
    # B (pxpP b) #
    if (!is.null(v.pa)) {B = qmat*Gamma}
    # kappa (px1 vector) #
    kappa = rep(1,p)

    s = 0

    while (s<S) {
        s = s+1
        if (s%%100==0) cat("no. of samples=",s,"\n")
        for (v in sample(1:p)) { # in random order
            # Update eta, A, kappa
            if (!is.null(v.pa)) {
                tempDat = dat.C - tcrossprodCpp(dat.P,B)
                no.de = rowSums(B!=0)
                bCb = sapply(1:p,function(v)sum(B[v,]^2/CC[v,]))
            } else{tempDat=dat.C
                no.de=bCb = rep(0,p)
            }
            up = updateUndirected(v=v,dat=tempDat,lambda=lambda,delta=delta,Alpha=A,eta=eta,kappa=kappa,pmat=pmat,no.de=no.de,bCb=bCb,no.tau=p)
            if (up$is.move) {
                eta=up$eta
                kappa = up$kappa
                A =up$Alpha
            }
            if (!is.null(v.pa)) {
                # Update Gamma, B, kappa
                v.ne = which(A[v,]!=0)
                v.ne.l = length(v.ne)
                vv.ne = c(v,v.ne)
                l.cl = length(vv.ne)
                if (v.ne.l>0) {
                    alpha = A[v,v.ne,drop=F]
                    y = c(dat.C[,v] - prod(dat.C[,v.ne,drop=F],alpha))
                    X = sapply(alpha,function(k) -k*dat.P,simplify=F)
                    X = cbind(dat.P,do.call(cbind,X))
                }else{
                    y = dat.C[,v]
                    X = dat.P
                    alpha=0
                }
                CCinv = (1/CC[vv.ne,,drop=F]) * kappa[vv.ne]
                CCinv = as.matrix(bdiag(sapply(1:l.cl,function(x)diag(CCinv[x,]),simplify=F)))
                tempB = B[vv.ne,,drop=F]
                tempGamma = Gamma[vv.ne,,drop=F]
                tempqmat = qmat[vv.ne,,drop=F]
                up = updateDirected(y=y,X=X,lambda=lambda,delta=delta,CCinv=CCinv,B=tempB,Gamma=tempGamma,kappa=kappa[v],qmat=tempqmat,v.no.ue =v.ne.l,v.aa=sum(alpha^2),no.tau=p)
                if (up$is.move){
                    Gamma[vv.ne,] = up$Gamma
                    kappa[v] = up$kappa
                    B[vv.ne,] = up$B
                }
            } #(!is.null(v.pa))
        }# for (v in sample(1:p))
        ### Store values
        if (s>burnin.S) {
            ss = s - burnin.S
            if (!is.null(v.pa)) {Gamma.list[,,ss] = Gamma; B.list[,,ss]=B}
            eta.list[,,ss] = eta
            A.list[,,ss] = A
            kappa.list[ss,] = kappa
        }
    }#while (s<S)
    return(list(Gamma=Gamma.list,eta=eta.list,A=A.list,B=B.list,kappa=kappa.list))
}


v.updateUndirected <- function(y,X,lambda,delta,Alpha,eta,kappa,pmat,no.de,bCb,no.tau) {
    # sampling alpha and kappa given eta
    # Input
    #   - y : response
    #   - X: design matrix
    #   - Afull: chol2inv(chol((crossprod(dat,tdat)+diag(p)*(1/lambda))))
    #   - lambda, delta: hyperparameters for kappa
    #   - eta: px1 eta matrix for v
    #   - kappa: px1 kappa vector (inverse covariance for all nodes)
    #   - pmat: pxp matrix for hyper parameters for eta
    #   - no.de : number of directed edges connected to each node (for updating kappa) (vector)
    #   - bCb : t(b) %*% inv(C) %*% b for updating kappa for each node (vector)
    #   - no.tau:  number of nodes in the current layer

    n = nrow(X)
    p = ncol(X)
    ### Graph movement ###
    if (sum(eta)==0) {is.swap=0
    }else{is.swap = rbinom(1,1,1/2)}
    is.move = FALSE
    if (!is.swap == 1) { #### ADD-DELETE
        w.ad = sample(1:p,1)
        new.eta = eta
        new.eta[w.ad] = 1-eta[w.ad]
        is.move = TRUE
    }else{ #### SWAP
        cand0 = which(eta==0)
        cand1 = which(eta==1)
        if (length(cand0)>0&length(cand1)>0) {
            w.ad1 = cand0[sample.int(length(cand0))[1]]
            w.ad2 = cand1[sample.int(length(cand1))[1]]
            new.eta = eta
            new.eta[w.ad1] = 1
            new.eta[w.ad2] = 0
            is.move = TRUE
        }
    }

    ### Sampling parameters ###
    if (is.move) {

        Afull = chol2inv(chol((crossprodCpp(X,X)+diag(p)*(1/lambda))))
        new.addr = which(new.eta==1)
        addr = which(eta==1)
        new.addr.inv =  which(new.eta==0)
        addr.inv = which(eta==0)
        A.new = giveA(Afull,which=new.addr.inv)
        A = giveA(Afull,which=addr.inv)

        #sample eta

        lhr1 = dmvnrm_arma(x = matrix(y,ncol=n),mean=rep(0,n),sigma = giveCov(X[,new.addr,drop=F],rep(lambda,length(new.addr)),kappa),log=TRUE)
        hyper1 = sum(new.eta*log(pmat) + (1-new.eta)*log(1-pmat))

        lhr2 = dmvnrm_arma(x = matrix(y,ncol=n),mean=rep(0,n),sigma =giveCov(X[,addr,drop=F],rep(lambda,length(addr)),kappa) ,log=TRUE)
        hyper2 = sum(eta*log(pmat) + (1-eta)*log(1-pmat))
        lhr = lhr1+hyper1-lhr2-hyper2


        if (log(runif(1))<lhr) {eta <- new.eta; A<-A.new}


        alpha = rep(0,p)
        resid = y
        if (!is.null(A)){
            # sample alpha
            addr = which(eta==1)
            if (length(addr)>0) {
                alpha[addr] = rmvnrm_arma(1,prodCpp(A,crossprodCpp(X[,addr,drop=F],as.matrix(y))),A/kappa)
            }
            # sample kappa
            resid = y -prodCpp(X,as.matrix(alpha))
        }
        kappa = rgamma(1,shape=(n+delta+no.tau-1+no.de+sum(alpha!=0))/2,rate=(lambda +sum(resid^2)+bCb+(lambda +no.tau-1)*sum(alpha^2) )/2)

    } #if (is.move)
    return(list(eta=eta,kappa=kappa,Alpha=alpha,is.move=is.move))
}


v.chaingraph <- function(v,chlist,palist,Y,eta.prob=0.1,gamma.prob=0.1,lambda,delta,burnin.S,inf.S) {
    # This function performs node-wise BANS regression model
    # Input
    # - v : target indice for node-wise regression
    # - chlist : list for components
    # - palist : list for parent set for each component
    # - Y : nxp data matrix

    t = which(unlist(lapply(chlist, function(x) v %in% x)))
    v.ch = chlist[[t]]
    stopifnot(length(v.ch)>2)
    v.pa = palist[[t]]
    S = burnin.S + inf.S
    dat.C = scale(Y[,v.ch,drop=F])
    n = nrow(dat.C)
    p = ncol(dat.C)
    v.addr = match(v,v.ch)
    pmat = rep(0.2,p-1)

    if (is.null(v.pa)) {
        Gamma.list = NULL
        B.list=NULL
    }else{
        dat.P = scale(Y[,v.pa,drop=F])
        pP = ncol(dat.P)
        Gamma.list = B.list = array(0,dim=c(p,pP,inf.S))
        CC = matrix(1/lambda,ncol=pP,nrow=p)# hyper parameter for b for nonzero gamma
        qmat = matrix(0.1,ncol(dat.C),ncol(dat.P))
    }
    eta.list =A.list= matrix(0,nrow=inf.S,ncol=p)
    kappa.list= rep(0,inf.S)

    ###########################################
    ####### Initialize all parameters #########
    ###########################################
    # eta (indicators for alpha) #
    eta = rep(0,p)
    eta[-v.addr] = rbinom(p-1,size=1,prob=eta.prob)
    # Gamma (indicators for b)#
    if (!is.null(v.pa)) {Gamma = matrix(rbinom(p*pP,size=1,prob=gamma.prob),nrow=p,ncol=pP)}
    # A (pxp alpha)$
    A = rep(0,p)
    A[-v.addr] = pmat*eta[-v.addr]
    # B (pxpP b) #
    if (!is.null(v.pa)) {B = qmat*Gamma}
    # kappa (px1 vector) #
    kappa = 1

    s = 0

    while (s<S) {
        s = s+1
        if (s%%100==0) cat("no. of samples=",s,"\n")
        # Update eta, A, kappa
        if (!is.null(v.pa)) {
            tempDat = dat.C - tcrossprodCpp(dat.P,B)
            no.de = sum(B[v.addr,]!=0)
            bCb = sum(B[v.addr,]^2/CC[v.addr,])
        } else{tempDat=dat.C
            no.de=bCb = rep(0,p)
        }
        up = v.updateUndirected(y=tempDat[,v.addr],X=tempDat[,-v.addr],lambda=lambda,delta=delta,Alpha=A[-v.addr],eta=eta[-v.addr],kappa=kappa,pmat=pmat,no.de=no.de,bCb=bCb,no.tau=p)
        if (up$is.move) {
            eta[-v.addr]=up$eta
            kappa = up$kappa
            A[-v.addr] =up$Alpha
        }
        if (!is.null(v.pa)) {
            # Update Gamma, B, kappa
            v.ne = which(A!=0)
            v.ne.l = length(v.ne)
            vv.ne = c(v.addr,v.ne)
            l.cl = length(vv.ne)
            if (v.ne.l>0) {
                alpha = A[v.ne,drop=F]
                y = c(dat.C[,v.addr] - prod(dat.C[,v.ne,drop=F],alpha))
                X = sapply(alpha,function(k) -k*dat.P,simplify=F)
                X = cbind(dat.P,do.call(cbind,X))
            }else{
                y = dat.C[,v.addr]
                X = dat.P
                alpha=0
            }
            CCinv = (1/CC[vv.ne,,drop=F]) * kappa
            CCinv = as.matrix(bdiag(sapply(1:l.cl,function(x)diag(CCinv[x,]),simplify=F)))
            tempB = B[vv.ne,,drop=F]
            tempGamma = Gamma[vv.ne,,drop=F]
            tempqmat = qmat[vv.ne,,drop=F]
            up = updateDirected(y=y,X=X,lambda=lambda,delta=delta,CCinv=CCinv,B=tempB,Gamma=tempGamma,kappa=kappa,qmat=tempqmat,v.no.ue =v.ne.l,v.aa=sum(alpha^2),no.tau=p)
            if (up$is.move){
                Gamma[vv.ne,] = up$Gamma
                kappa= up$kappa
                B[vv.ne,] = up$B
            }
        } #(!is.null(v.pa))
        ### Store values
        if (s>burnin.S) {
            ss = s - burnin.S
            if (!is.null(v.pa)) {Gamma.list[,,ss] = Gamma; B.list[,,ss]=B}
            eta.list[ss,] = eta
            A.list[ss,] = A
            kappa.list[ss] = kappa
        }
    }#while (s<S)
    return(list(Gamma=Gamma.list,eta=eta.list,A=A.list,B=B.list,kappa=kappa.list))
}


ch.chaingraph.str<- function(v.ch,v.pa,Y,G,lambda,delta,burnin.S,inf.S) {
# This function fits regression model for a chain component (structured estimation)
# Input
# - v.ch : indices for the target chain component
# - v.pa : indices for parents set
# - Y : nxp data matrix
# - G : structure
	stopifnot(length(v.ch)>2)
	S = burnin.S + inf.S
	dat.C = scale(Y[,v.ch,drop=F])
	n = nrow(dat.C)
	p = ncol(dat.C)

	if (is.null(v.pa)) {
		B.list=NULL
	}else{
		dat.P = scale(Y[,v.pa,drop=F])
		pP = ncol(dat.P)
		B.list = array(0,dim=c(p,pP,inf.S))
		CC = matrix(1/lambda,ncol=pP,nrow=p)# hyper parameter for b for nonzero gamma
	}
    A.list= array(0,dim=c(p,p,inf.S))
	kappa.list = matrix(0,nrow=inf.S,ncol=p)
###########################################
####### Initialize all parameters #########
###########################################
# eta (indicators for alpha) #
	eta = G[v.ch,v.ch]
# Gamma (indicators for b)#
	if (!is.null(v.pa)) {Gamma = G[v.ch,v.pa]}
# A (pxp alpha)$
	A = 0.2*eta
# B (pxpP b) #
	if (!is.null(v.pa)) {B = 0.2*Gamma}
# kappa (px1 vector) #
	kappa = rep(1,p)

	s = 0

    while (s<S) {
        s = s+1
        if (s%%1000==0) cat("no. of samples=",s,"\n")
        for (v in sample(1:p)) { # in random order
# Update eta, A, kappa
        if (!is.null(v.pa)) {
                tempDat = dat.C - tcrossprodCpp(dat.P,B)
                    no.de = rowSums(B!=0)
                    bCb = sapply(1:p,function(v)sum(B[v,]^2/CC[v,]))
        } else{tempDat=dat.C
            no.de=bCb = rep(0,p)
        }
            up = updateUndirected.str(v=v,dat=tempDat,lambda=lambda,delta=delta,Alpha=A,eta=eta,kappa=kappa,no.de=no.de,bCb=bCb,no.tau=p)
			kappa = up$kappa
			A =up$Alpha

            if (!is.null(v.pa)) {
# Update Gamma, B, kappa
                v.ne = which(A[v,]!=0)
                v.ne.l = length(v.ne)
                vv.ne = c(v,v.ne)
                l.cl = length(vv.ne)
                if (length(v.ne)>0) {
                    alpha = A[v,v.ne,drop=F]
                    y = c(dat.C[,v] - prod(dat.C[,v.ne,drop=F],alpha))
                    X = sapply(alpha,function(k) -k*dat.P,simplify=F)
                    X = cbind(dat.P,do.call(cbind,X))
                }else{
                    y = dat.C[,v]
                    X = dat.P
                    alpha=0
                }
                CCinv = (1/CC[vv.ne,,drop=F]) * kappa[vv.ne]
                CCinv = as.matrix(bdiag(sapply(1:l.cl,function(x)diag(CCinv[x,]),simplify=F)))
                tempB = B[vv.ne,,drop=F]
                tempGamma = Gamma[vv.ne,,drop=F]
                up = updateDirected.str(y=y,X=X,lambda=lambda,delta=delta,CCinv=CCinv,B=tempB,Gamma=tempGamma,kappa=kappa[v],v.no.ue =v.ne.l,v.aa=sum(alpha^2),no.tau=p)
				kappa[v] = up$kappa
				B[vv.ne,] = up$B
            } #(!is.null(v.pa))
        }# for (v in sample(1:p))
### Store values
		if (s>burnin.S) {
			ss = s - burnin.S
			if (!is.null(v.pa)) {B.list[,,ss]=B}
			A.list[,,ss] = A
			kappa.list[ss,] = kappa
		}
    }#while (s<S)
	return(list(A=A.list,B=B.list,kappa=kappa.list))
}



updateUndirected.str <- function(v, dat,lambda,delta,Alpha,eta,kappa,no.de,bCb,no.tau) {
# sampling alpha and kappa given eta
# Input
#   - v: current node
#   - dat: nxp data matrix
#   - Afull: chol2inv(chol((crossprod(dat,tdat)+diag(p)*(1/lambda))))
#   - lambda, delta: hyperparameters for kappa
#   - eta: pxp eta matrix
#   - kappa: px1 kappa vector (inverse covariance for all nodes)
#   - pmat: pxp matrix for hyper parameters for eta
#   - no.de : number of directed edges connected to each node (for updating kappa) (vector)
#   - no.de : t(b) %*% inv(C) %*% b for updating kappa for each node (vector)
#   - no.tau:  number of nodes in the current layer

	addr = which(eta[v,]==1)
	y = dat[,v]
	n = length(y)
	p = length(addr)
	alpha = rep(0,ncol(dat))
	if (p>0) {
		X = dat[,addr,drop=F]
		A= chol2inv(chol((crossprodCpp(X,X)+diag(p)*(1/lambda))))
		alpha[addr] = rmvnrm_arma(1,prodCpp(A,crossprodCpp(X,as.matrix(y))),A/kappa[v])
	}
	Alpha[v,] = alpha
	resid = dat[,v] -prodCpp(dat,as.matrix(alpha))
	kappa[v] = rgamma(1,shape=(n+delta+no.tau-1+no.de[v]+sum(alpha!=0))/2,rate=(lambda +sum(resid^2)+bCb[v]+(lambda +no.tau-1)*sum(alpha^2) )/2)
    return(list(kappa=kappa,Alpha=Alpha))
}


updateDirected.str <- function(y,X,lambda,delta,CCinv,B,Gamma,kappa,v.no.ue,v.aa,no.tau) {
# sampling beta and kappa given gamma
# Input
#   - y: nx1 response vector corresponding to node v
#   - X: nx(p*pP) covariate vector including parents of v and a function of parents of C(v)
#   - lambda, delta: hyperparameters for kappa
#   - CCinv: pxpP hyper parameter matrix for B
#   - B : pxpP current B (permuted for response y)
#   - Gamma: pxpP current gamma (permuted for response y)
#   - kappa: a scalar for kappa for v (inverse covariance)
#   - v.no.ue: number of undirected edges for v (scalar) (for sampling kappa)
#   - v.aa: t(alpha) %*% alpha (scalar) (for sampling kappa)
#   - no.tau:  number of nodes in the current layer

    n = length(y)
    pP = ncol(B)
    p = nrow(B)
    vecpP = prod(dim(B))

### Sampling parameters ###
	vB = c(t(B))
	vGamma = c(t(Gamma))
	CC = 1/diag(CCinv)
	addr = which(vGamma==1)
	addr.inv = which(vGamma==0)

	if (pP*p>n) {
		XC =  X%*%diag(CC)
		tmp1 =chol2inv(chol(diag(n) + kappa * tcrossprodCpp(XC,X)))
		Afull = diag(CC)-kappa*crossprodCpp(XC,prodCpp(tmp1,XC))
	}else{Afull = chol2inv(chol(kappa*crossprodCpp(X,X)+CCinv))}
	A = giveA(Afull,which=addr.inv)

# sample B
	vB = rep(0,length(vGamma))
	resid = y
	if (!is.null(A)){
		ww = which(vGamma==1)
		if (length(ww)>0) {
			vB[ww]= rmvnrm_arma(1,kappa * prodCpp(A,crossprodCpp(X[,ww,drop=F],as.matrix(y))),A)
		}
		resid = y-prodCpp(X,as.matrix(vB))
	}
        kappa = rgamma(1,shape=(n+delta+no.tau-1+v.no.ue+sum(vB!=0))/2,rate=(lambda+sum(resid^2)+sum(vB^2/diag(CCinv))+(lambda +no.tau-1)*v.aa )/2)
    B = matrix(vB,ncol=pP,byrow=T)
    return(list(kappa=kappa,B=B))
}


