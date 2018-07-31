.loadPackage <- function(pkg){
  
  if(length(grep(pkg, search())) >0 ){
    return(invisible(TRUE))
  }
  
  result <- tryCatch(suppressPackageStartupMessages(
    require(pkg, quietly=TRUE, character.only=TRUE)
  ))
  if (!result){
    stop(paste("Please install package",pkg))
  }
  return(invisible(result))
}

mydist<-function(loc1,loc2){
X=as.matrix(loc1)
Y=as.matrix(loc2)
dist=outer(seq(length(X)),seq(length(Y)), Vectorize(function(i, j) abs(X[i]-Y[j])))
return(dist)
}

read2dense<-function(hicfile,N,resolution){
  # read hicmatrix
  if(length(hicfile)==1){
  triple_ = read.table(hicfile,stringsAsFactors =FALSE)
  triple = as.matrix(triple_)
  }else{
  triple= as.matrix(hicfile)
  }

 if (nrow(triple)==ncol(triple)){
    im=triple
    N=nrow(triple)
  }else{
    if(missing(N)){
     N=max(triple[,1])
     N=max(c(triple[,2],N))/resolution+1;
    }
    im = zeros(N,N)
    im[triple[,1]/resolution+1+N*(triple[,2]/resolution)]=triple[,3]
    im = im+t(im)
    im = im - diag(diag(im))
  
  }
  return(im)
}

KRnorm <-function(im){
	# find gap
   pos=which(colSums(im>0)==0)
   pos2=which(colSums(im>0)!=0)
   im2=im;
   im2=im2[-pos,]
   im2=im2[,-pos]
   n=nrow(im2)
   A=matrix(1,41,n)
   for (i in 1:20){
   A[21+i,1:(n-i)]=Diag(im2,-i)
   A[21-i,(i+1):n]=Diag(im2,i)
   }
   gapidx=which(colSums(A>0)<35)
   gapidx=c(pos,pos2[gapidx])
#         sumim=colSums(im>0)
#	cutoff = as.numeric(quantile(sumim[which(sumim>0)],0.05))
#	gapidx=which(colSums(im)<cutoff)


	if(sum(round(im)!=im)==0){
		    message('Your input matrix is raw matrix, perform KR normalization');
	    imnew_ = bnewt2(im)
	        imnew = as.matrix(imnew_$Anew)
	      }else{
		     message('Your input matrix is normalized matrix, skip KR normalization');
	          imnew=im
		    }

	return(list(imKR=imnew,gapidx=gapidx))
}


calcu_aRI <-function(imnew,gapidx,wd,wdsize){
  # wdsize=6 is the best in most exmaples we tes
	N=nrow(imnew)
	oridx=as.matrix(c(1:N))
	oridx = oridx[-gapidx]
	imnew = imnew[-gapidx,]
	imnew = imnew[,-gapidx]

	# calculate RIP
	tmpx = log(1+imnew)
	tmpn = nrow(tmpx)
	tads = zeros(wdsize,tmpn)

	for (w in( wd:(wd+wdsize-1) )){
		  xx = (w+1):(tmpn-w-1)
	  sumup = 0
	    for (i in (-w:-1))
		        for (j in ((i+1):0))
				      sumup = sumup+tmpx[xx+i+ tmpn*( xx+j-1)]

	    sumdown = 0
	      for (i in 1:w)
		          for (j in (i+1):(w+1))
				        sumdown = sumdown+tmpx[xx+i+ tmpn*( xx+j-1)]

	      sumrect = 0
	        for (i in (-w+1):0)
			    for (j in (1:(w+i)))
				          sumrect = sumrect+tmpx[xx+i+ tmpn*( xx+j-1)]

	        tadscores = (sumup+sumdown-sumrect)/(sumup+sumdown+sumrect)
		  tadscores = cbind (t(rep(tadscores[1],w)), t(tadscores) ,t(rep(tadscores[length(tadscores)],w+1)))
		  tads[w-wd+1,] = tadscores
	}


	# mean
	tadscores = colMeans(tads)

	#aRI
	tmp=NaN*zeros(N,1)
	tmp[oridx]=tadscores
	return(tmp)
}



getNorm <- function(replicates = list(),plt=FALSE,fout="M-A plot.png"){

    if(length(replicates) < 2)
    stop("At least two replicates should be provided")
  nbBins <- dim(replicates[[1]])
  smry1=Matrix::summary(replicates[[1]])
  smry2=Matrix::summary(replicates[[2]])
  mrg <- merge(smry1, smry2, by=c("i","j"), all=TRUE, suffixes = c(1,2))
  rm(smry1,smry2)

  mrg$x1[is.na(mrg$x1)] <- 0
  mrg$x2[is.na(mrg$x2)] <- 0
  mrg$x1 <- mrg$x1+1
  mrg$x2 <- mrg$x2+1

  M <- log2( mrg$x1 / mrg$x2)
  A <- 0.5 * log2(mrg$x1 * mrg$x2)

  o <- order(A)
  A <- A[o]; M <- M[o]

  fit <- loessFit(M,A, iterations = 5)

  if(plt){
    message(paste("Plotting M-A plot into file",fout))
    png(fout, width=1024, height=1024)
    plot(A,M, main="M-A plot between replicates")
    lines(A,fit$fitted, col='red',lwd=2)
    dev.off()
  }


  mrg <- mrg[o,]

  replicates[[1]] <- sparseMatrix(i= mrg$i, j=mrg$j, x= 2^(log2(mrg$x1) - 0.5 * fit$fitted),dims=nbBins)
  replicates[[2]] <- sparseMatrix(i= mrg$i, j=mrg$j, x= 2^(log2(mrg$x2) + 0.5 * fit$fitted),dims=nbBins)
  return(replicates)
}



chrpeaks <-function(hicfile,resolution,N,mind,wd,wdsize,chr)
{
  # Detect local maximum peaks on a single chromosome
  # Implemented by
  # Fengling Chen,Guipeng Li
  # Tsinghua University
  # https://github.com/ChenFengling

  # wdsize=6 is the best in most exmaples we test
  chrsize=N
  N=ceiling(N/resolution)
  im=read2dense(hicfile,N,resolution)
  out=KRnorm(im)
  imnew=out$imKR
  gapidx=out$gapidx

message('find local maximum')
  oridx=as.matrix(c(1:N))
  oridx = oridx[-gapidx]
  imnew = imnew[-gapidx,]
  imnew = imnew[,-gapidx]

  # calculate RIP
  tmpx = log(1+imnew)
  tmpn = nrow(tmpx)
  tads = zeros(wdsize,tmpn)

  for (w in( wd:(wd+wdsize-1) )){
    xx = (w+1):(tmpn-w-1)
    sumup = 0
    for (i in (-w:-1))
      for (j in ((i+1):0))
        sumup = sumup+tmpx[xx+i+ tmpn*( xx+j-1)]

    sumdown = 0
    for (i in 1:w)
      for (j in (i+1):(w+1))
        sumdown = sumdown+tmpx[xx+i+ tmpn*( xx+j-1)]

    sumrect = 0
    for (i in (-w+1):0)
      for (j in (1:(w+i)))
        sumrect = sumrect+tmpx[xx+i+ tmpn*( xx+j-1)]

    tadscores = (sumup+sumdown-sumrect)/(sumup+sumdown+sumrect)
    tadscores = cbind (t(rep(tadscores[1],w)), t(tadscores) ,t(rep(tadscores[length(tadscores)],w+1)))
    tads[w-wd+1,] = tadscores
  }


  # mean
  tadscores = colMeans(tads)
 if(length(which(is.nan(tadscores)))!=0){
         pos=which(is.nan(tadscores))
	 tadscores=tadscores[-pos]
	 oridx=oridx[-pos]
         gapidx=sort(c(gapidx,oridx[pos]))
           tmpn=length(oridx)
 }

  # calculate envelope
  temp = envelope(tadscores)
  ylower = as.matrix(temp[4])

  temp = envelope(ylower)
  ylower2 = as.matrix(temp[4])

  peakscore=tadscores-ylower2
  LRI=peakscore
  #find peaks

  temp = findpeaks(tadscores,minpeakdistance=mind)
  TAD_boundaries = temp[,2]

  #rescue  neighbouring peak with  high score
  temp = findpeaks(tadscores)
  TAD_boundaries2 = temp[,2]

  highcut=sort(peakscore[TAD_boundaries],decreasing=T)

  highcut=highcut[round(tmpn/(500000/resolution))]

  logical1=peakscore>highcut
  logical2=zeros(1,tmpn)
  logical2[TAD_boundaries2]=1
  logical3=zeros(1,tmpn)
  logical3[TAD_boundaries]=1
  logical1=t(as.matrix(logical1))
  rescue=which(logical1==1&logical2==1&logical3==0)
  TAD_boundaries=cbind(t(TAD_boundaries),rescue)
  peakscore=peakscore[TAD_boundaries]
  peakavgRI=tadscores[TAD_boundaries]
  boundary=rbind(TAD_boundaries,peakscore,peakavgRI)

  # change to original  boundary
  boundary[1,]=oridx[boundary[1,]]

  # remove repeat region near
  temp = distmat(t(t(gapidx)),t(t(boundary[1,])))
  rRa = apply(temp,2,min)
  rpflt = which(rRa>wd)
  boundary=boundary[,rpflt]
  peakidx=order(boundary[2,],decreasing = T)
  boundary=boundary[,peakidx]
  boundary=boundary[,!is.na(boundary[2,])]
  boundary=rbind(rep(chr,1,ncol(boundary)),boundary[1,]*resolution-resolution/2,boundary[1,]*resolution+resolution/2,boundary[2:nrow(boundary),])
  boundary=t(boundary)

  #aRI
  tmp=NaN*zeros(nrow(im),1)
  tmp[oridx]=tadscores
  tadscores=tmp
  Ns=nrow(im)-1
  allscore=rbind(rep(chr,1,Ns),(1:Ns)*resolution-resolution/2,c((1:(Ns-1))*resolution+resolution/2,chrsize),t(tadscores[1:Ns,1]))
  allscore=t(allscore)

  #LRI
  tmp=NaN*zeros(nrow(im),1)
  tmp[oridx]=LRI;
  LRI=tmp
  Ns=nrow(im)-1
  allLRI=rbind(rep(chr,1,Ns),(1:Ns)*resolution-resolution/2,c((1:(Ns-1))*resolution+resolution/2,chrsize),t(LRI[1:Ns,1]))
  allLRI=t(allLRI)
  return(list(LRI = allLRI,aRI =allscore ,localmax=boundary))
}

chrpeaks_diff <-function(hicfile,resolution,N,mind,wd,wdsize,chr)
{
	chrsize=N
	N=ceiling(N/resolution)
	##  step.1 calucate aRI score
        n1=length(hicfile[[1]])
        n2=length(hicfile[[2]])
	sampleA=list(NULL)
	sampleB=list(NULL)
	gapidxA=list(NULL)
	gapidxB=list(NULL)
	sumKR=c()
	aRI=matrix(0,N,n1+n2)
	mergedA=matrix(0,N,N)
	mergedB=matrix(0,N,N)
	## KR norm + depthnorm
	# sample A
	
	message("step.1 calucate aRI score on each replicate")
	for (i in 1:n1){
	    out=read2dense(paste(hicfile[[1]][i],'/chr',chr,'.matrix',sep=""),N,resolution)
	    mergedA=mergedA+out
	    out=KRnorm(out)
	    sumKR[i]=sum(out$imKR)
	    sampleA[[i]]=Matrix(out$imKR/sumKR[i]*sumKR[1],sparse =TRUE)
	    gapidxA[[i]]=out$gapidx
	}
	
	sampleA=getNorm(sampleA,plt=1,fout ="sampleA_MA_before.png" )

	for (i in 1:n1){
		  aRI[,i]=calcu_aRI(as.matrix(sampleA[[i]]),gapidxA[[i]],wd,wdsize)
	}

	# sample B
	for (i in 1:n2){
      out=read2dense(paste(hicfile[[2]][i],'/chr',chr,'.matrix',sep=""),N,resolution)
	    mergedB=mergedB+out
	    out=KRnorm(out)
	    sumKR[n1+i]=sum(out$imKR)
	    sampleB[[i]]=Matrix(out$imKR/sumKR[n1+i]*sumKR[1],sparse=TRUE)
	    gapidxB[[i]]=out$gapidx
		     }
	sampleB=getNorm(sampleB,plt=1,fout ="sampleB_MA_before.png" )
	for (i in 1:n2){
		  aRI[,n1+i]=calcu_aRI(as.matrix(sampleB[[i]]),gapidxB[[i]],wd,wdsize)
	}
        message("step2. calculate CDBs in merged samples")
	## step2.calculate CDBs in merged samples
	resultA=chrpeaks(mergedA,resolution,chrsize,mind,wd,wdsize,chr)
	resultB=chrpeaks(mergedB,resolution,chrsize,mind,wd,wdsize,chr)
	return(list(aRI=aRI,resultA=resultA,resultB=resultB))
}


EScutoff <-function(maxall,motiffile,outdir){
# Decide a boundary cut off by  CTCF notif enrichment
#
# boundary: boundary with  LRI score
# motiffile: CTCF sites
# resolution: resolution of HiC matrix
#
# Implemented by  Fengling Chen
# Tsinghua University
# https://github.com/ChenFengling/
  motiffile_ = read.table(motiffile)
  motiffile = as.matrix(motiffile_)

  peakidx =order(maxall[,4],decreasing = T)
  maxall=maxall[peakidx,]
  chr=unique(maxall[,1])
  allpos=zeros(nrow(maxall),1)
  for (i in 1:length(chr)){
    motif=motiffile[which(motiffile[,1]==chr[i]),]
    chrpos=which(maxall[,1]==chr[i])
    localmax=maxall[chrpos,]
    X=as.matrix(motif[,2])
    Y=as.matrix((localmax[,2]+localmax[,3])/2)
    dist=outer(seq(length(X)),seq(length(Y)), Vectorize(function(i, j) abs(X[i]-Y[j])))
    pos = apply(dist,2,min) <20000
    allpos[chrpos]=pos
  }
  pcoe=maxall[which(allpos==1),4]/sum(maxall[which(allpos==1),4])
  ncoe=1/sum(allpos==0)
  allpos[which(allpos==0)]= -ncoe
  allpos[which(allpos==1)]= pcoe
  score= cumsum(allpos)
  I=which.max(score)
  boundary=maxall[1:I,]
  return(list(boundary=boundary,score=score))
}



envelope <-function(sig, method='linear')
{
# Find upper and lower envelopes of a given signal
# The idea is from Envelope1.1 by Lei Wang, but here it works well when the signal contains
# successive equal samples and also includes first and last samples of the signal in the envelopes.
# inputs:
#   sig: vector of input signal
#   method: method of interpolation (defined as in interp1)

# outputs:
  #   upperenv: upper envelope of the input signal
#   lowerenv: lower envelope of the input signal

  upperind = which(diff(sign(diff(sig))) < 0) + 1
  lowerind = which(diff(sign(diff(sig))) > 0) + 1
  f = 1
  l = length(sig)


upperind = rbind(f,as.matrix(upperind),l)
lowerind = rbind(f,as.matrix(lowerind),l)

xi = f : l
upperenv = approx(upperind, sig[upperind], xi, method)
lowerenv = approx(lowerind, sig[lowerind], xi, method)
return(data.frame(upperenv = upperenv,lowerenv= lowerenv))
}



bnewt2= function(A0,tol=1e-5,delta=0.05,Delta=2,fl=0)
{
n0 = nrow(A0)
KR0 = which(colSums(A0)>100)
A = A0[KR0,KR0]

# BNEWT A balancing algorithm for symmetric matrices
#
#X = BNEWT(A) attempts to find a vector X such that
#diag(X)*A*diag(X) is close to doubly stochastic. A must
#be symmetric and nonnegative.
#
#X0: initial guess. TOL: error tolerance.
#delta/Delta: how close/far balancing vectors can get
#to/from the edge of the positive cone.
#We use a relative measure on the size of elements.
#FL: intermediate convergence statistics on/off.
#RES: residual error, measured by norm(diag(x)*A*x - e).
#Initialise
n = nrow(A)
e = ones(n,1)
x0 = e
res = NaN
# Inner stopping criterion parameters.
g=0.9
etamax = 0.05
eta = etamax
stop_tol = tol*.5

x = x0; rt = tol^2; v = x*(A%*%x); rk = 1 - v;

rho_km1 = t(rk)%*%rk; rout = rho_km1; rold = rout;
MVP = 0;
i = 0; # Outer iteration count.
if (fl == 1) {fprintf('it in. it res\n')}
while (rout > rt) { # Outer iteration
  i = i + 1; k = 0; y = e;
  innertol = max(eta^2*rout,rt)
  while (rho_km1 > innertol){ #Inner iteration by CG
    k = k + 1
    if (k == 1){
      Z = rk/v; p=Z; rho_km1 = t(rk)%*%Z;
    }else{
      beta=rho_km1/rho_km2
      p=Z + as.numeric(beta)*p
    }
    # Update search direction efficiently.
    w = x*(A%*%(x*p)) + v*p
    alpha = rho_km1/(t(p)%*%w)
    ap = as.numeric(alpha)*p
    # Test distance to boundary of cone.
    ynew = y + ap;
    if (min(ynew) <= delta){
      if (delta == 0) {break}
      ind = which(ap < 0)
      gamma = min((delta - y[ind])/ap[ind])
      y = y + gamma*ap;
      break
    }
    if (max(ynew) >= Delta){
      ind = which(ynew > Delta);
      gamma = min((Delta - y[ind])/ap[ind])
      y = y + gamma*ap
      break
    }
    y = ynew
    rk = rk - as.numeric(alpha)*w; rho_km2 = rho_km1;
    Z = rk/v; rho_km1 = t(rk)%*%Z;
  }
    x = x*y; v = x*(A%*%x);
    rk = 1 - v; rho_km1 = t(rk)%*%rk; rout = rho_km1;
    MVP = MVP + k + 1;

    # Update inner iteration stopping criterion.
    rat = rout/rold; rold = rout; res_norm = sqrt(rout);
    eta_o = eta; eta = g%*%rat;
    if (g%*%eta_o^2 > 0.1){
    eta = max(eta,g*eta_o^2)
    }

    eta = max(min(eta,etamax),stop_tol/res_norm);
    if (fl == 1){
      fprintf('%3d %6d %.3e \n',i,k, r_norm);
      res=r_norm
    }
}
    KRnorm = ones(n0,1)
    KRnorm[KR0] = x
    s = sum(A0)/n0
    KRnorm = KRnorm%*%sqrt(s)
    Anew = A0*(KRnorm%*%t(KRnorm))
    return(list(Anew=Anew,KRnorm=KRnorm,res=res))
}
