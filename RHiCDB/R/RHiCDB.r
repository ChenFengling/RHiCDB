RHiCDB=function(hicfile,resolution,chrsizes,ref='no',outdir,mind,wd,wdsize){
  #'@title Detect CDBs and differential CDBs on Hi-C heatmap.
  #'   
  #'@description HiCDB using Hi-C  contact matrix to detect Hi-C contact domain boundaries(CDBs).It outputs annotated CDBs, differential CDBs  on the chosen options
  #'
  #'@details 
  #'    A. Possible outputs 
  #'     
  #'    1.CDB.txt
  #'    
  #'    2.localmax.txt: all the local maximum peaks detected before cutoff
  #'    decision. User can decide custum CDB cutoff upon this file.
  #'    
  #'    3.EScurve.png: CTCF motif enrichment on ranked local maximum peaks.
  #'    
  #'    4.aRI.txt: average RI score for each  genomic bin.  
  #'    
  #'    5.LRI.txt: LRI score for each genomic bin.
  #'    
  #'     B. default value for 'mind','wd' on different resolution
  #'   
  #'     resolution mind  wd  wdsize
  #'    
  #'    10k         4    3    6
  #'    
  #'    40k         2    1    6
  #'    
  #'    5k          8    6    8      
  #'    
  #'    C. HiCDB will perform a KR normlization if the data is raw counts.
  #'@param hicfile hicfile is the directory of the intra-chromosome Hi-C matrixes with 
  #'    sparse or dense format. The intra-chromosome matrix must be named as "chr+number.matrix" 
  #'    according to the chromosome order like 'chr1.matrix','chr2.matrix',..., 'chr23.matrix'. 
  #'    As HiCDB matches "chr*.matrix" to recognize the Hi-C matrix, avoid to use the "chr*.matrix" as 
  #'    the name of other files. The intra-chromosome matrix could be in a dense (a NxN matrix) or sparse (a Kx3 table,Rao et al.) format. 
  #'
  #'    If you want to detect CDB on one sample,set hicfile as 'SAMPLE_DIR'. 
  #'    If ref is not set, this function will output all the local maximum peaks. 
  #'    If ref is set, this function will output local maximum peaks and final CDBs.
  #'    
  #'    If you want to detect differential CDBs, ref is required to decide the cut off on CDB detection.
  #'    If you don't have replicate, set hicfile as list('SAMPLE1','SAMPLE2'). This function will first perform CDB detection on each sample and then compare the difference between their final CDBs by intersection. If replicates is provided, set hicfile as list(c('SAMPLE1_rep1','SAMPLE1_rep2'),c('SAMPLE2_rep1','SAMPLE2_rep2')). The function will find CDBs on each sample with  merged Hi-C matrix, calculate aRI score on each replicates, then decide a CDB as differetial or not by statistical test on  aRI scores of each CDB. 
  #'    
  #'    If ref is 'hg38' or 'hg19', CDBs will also be annotated as conserved or not conserved.
  #'    
  #'    
  #' @param resolution  resolution of Hi-C matrix. This is required.
  #'
  #' @param chrsizes  Ordered chromosome sizes of the genome. Optional setting is 'hg19', 'hg38', 'mm9', 'mm10' 
  #'       or any other chromosome size files which can be generated following the instructions in annotation/README.md. 
  #'       This is required.
  #'       
  #'@param ref   reference CTCF motif locs on the genome. If it is set, the output will use
  #'    the GSEA-like methods to decide the cutoff. Default is 'no'. Choices are :
  #'       'no'
  #'       'hg19'
  #'       'hg38'
  #'       'mm9'
  #'       'mm10'
  #'        or other customfile  for example 'genome.txt' made from utility/motifanno.sh
  #'         Example for 'genome.txt':
  #'             #'chr   motifcenterlocus
  #'             10	15100928
  #'             10	15188593
  #'
  #' @param outdir   The output directory. Default will be the directory of the first sample.
  #'
  #' @param mind    Minimum local maximum peak distance (measured by bin), or
  #'                 minimum separation between local maximum peaks, specified
  #'                as a positive integer scalar. Use this argument to have findpeaks
  #'                ignore small peaks that occur in the neighborhood of a larger peak.
  #'
  #'@param wd       The smallest window sizes.
  #'
  #'@param wdsize   The number of different window size. The whole window size
  #'                 scale will be wd:(wd+wdsize).Default will be 6.
  #'    
  #'@author
  #'   Implemented by Fengling Chen
  #' 
  #'   Any suggestions and remarks might be addressed to Fengling Chen:cfl15@mails.tsinghua.edu.cn
  #'@export
  #'@examples
  #'      1. Output all the local maximum peaks and let customers to decide the cutoff.
  #'      HiCDB('sample1/',10000,chrsizes='custom_chrsizes.txt');
  #'      HiCDB('sample1/',10000,chrsizes='custom_chrsizes.txt',outdir='sample1/outputs/');
  #'      2. Use GSEA-like methods to decide the cutoff
  #'      HiCDB('sample1/',10000,chrsizes='hg19',ref='hg19');
  #'      HiCDB('sample1/',10000,chrsizes='custom_chrsizes.txt',ref='custom_motiflocs.txt');
  #'      3. To detect differential CDBs
 #'       HiCDB(list('sample1','sample2'),10000,'hg19',ref='hg19');
 #'       HiCDB(list(c("sample1_rep1","sample1_rep2"),c("sample2_rep1","sample2_rep2")), 
  #'	 + 10000,'hg19',ref='hg19');


  
  options(scipen = 20)
  .loadPackage("pracma") 
  .loadPackage("Matrix") 
  .loadPackage("limma") 
  ################################################################################################################
  
  ##############################  check input ######################################
  
  ##################################################################################################################
if (length(hicfile)==1){
    print("One sample is provided, detect CDBs!")
    FLAG=1
  }else if(length(hicfile)==2){
    if(ref=='no'){
      stop("You should set ref in differential CDB detection!")
    }
    if(length(hicfile[[1]])==1){
      FLAG=2
      print('differential CDB called by intersection!')
    }else{
      FLAG=3
      print('differential CDB called by statistical test!')
    }
  }else{
    stop("You should provided samples in a list with length 1 or 2!")
  }
  
if (is.null(resolution)){
    stop("resolution is required!")
}

if(!is.null(chrsizes)){  
if(chrsizes=='hg38'){
  chrsizes=system.file("extdata", "hg38.chrom.sizes.txt", package="RHiCDB")
  }else if(chrsizes=='hg19'){
    chrsizes=system.file("extdata", "hg19.chrom.sizes.txt", package="RHiCDB")
  }else if(chrsizes=='mm9'){
    chrsizes=system.file("extdata", "mm9.chrom.sizes.txt", package="RHiCDB")
  }else if(chrsizes=="mm10"){
    chrsizes=system.file("extdata", "mm10.chrom.sizes.txt", package="RHiCDB")
  }else{
    chrsizes=chrsizes
  }
}else{
  stop("chrsizes is required!")
}
chrsizes=read.table(chrsizes)
chrsizes=chrsizes[,1]
#chrsizes=floor(chrsizes[,1]/resolution)+1

if (missing(outdir)){
  outdir=hicfile
}else if(!dir.exists(outdir)){
       dir.create(outdir,recursive=TRUE)
    }

 
if(missing(wd)){
    if(resolution==40000){
          wd=1
    }else if(resolution==10000){
          wd=3
    }else if(resolution==5000){
          wd=6
    }else{
      stop("please set parameter wd under this resolution")
    }
  }
  
  if(missing(wdsize)){
    if(resolution==40000){
      wdsize=6
    }else if(resolution==10000){
      wdsize=6
    }else if(resolution==5000){
      wdsize=8
    }else{
      stop("please set parameter wdsize under this resolution")
    }
  }
  
if(missing(mind)){
    mind=40000/resolution
  }
  
if(ref!='no'){
  if(ref=='hg38'){
    motiffile=system.file("extdata", "CTCF_hg38_site.txt", package="RHiCDB")
  }else if(ref=='hg19'){
    motiffile=system.file("extdata", "CTCF_hg19_site.txt", package="RHiCDB")
  }else if(ref=='mm9'){
    motiffile=system.file("extdata", "CTCF_mm9_site.txt", package="RHiCDB")
  }else if(ref=="mm10"){
    motiffile=system.file("extdata", "CTCF_mm10_site.txt", package="RHiCDB")
  }else{
    motiffile=ref
  }
}
################################################################################################################

##############################  CDB detection ######################################

##################################################################################################################

if(FLAG==1){
  tmp=list.files(path=hicfile,pattern='chr.*\\.matrix')
  totalchr=length(tmp)
  allpeaks=as.data.frame(matrix(numeric(0),ncol=4))
  allLRI=as.data.frame(matrix(numeric(0),ncol=4))
  allaRI=as.data.frame(matrix(numeric(0),ncol=4))
  print('* Detect all local maximum peaks on each chromosome.')
  
  for( i in 1:totalchr){
    print(paste('Processing chr',i,sep=""))
    hicmap=paste(hicfile,'/chr',i,'.matrix',sep="")
    out=chrpeaks(hicmap,resolution,chrsizes[i],mind,wd,wdsize,i)
    allpeaks=rbind(allpeaks,out$localmax)
    allLRI=rbind(allLRI,out$LRI)
    allaRI=rbind(allaRI,out$aRI)
  }
  peakidx=order(allpeaks[,4],decreasing = TRUE)
  allpeaks=allpeaks[peakidx,]
  print(' ')
  print('* Write outputs: localmax.txt.')
  write.table(allpeaks,paste(outdir,'/localmax.txt',sep=""),sep="\t",col.names=F,row.names=F,quote=F)
  write.table(allaRI,paste(outdir,'/aRI.txt',sep=""),sep="\t",col.names=F,row.names=F,quote=F)
  write.table(allLRI,paste(outdir,'/LRI.txt',sep=""),sep="\t",col.names=F,row.names=F,quote=F)
  if(ref!="no"){
  print(' ')
  print('* Decide cut off using enrichmentscore. Write outputs:EScurve.png.')
  out=EScutoff(allpeaks,motiffile,outdir)
  CDB=out$boundary
  score=out$score
  pdf(paste(outdir,'/EScurve.pdf',sep=""))
  plot(score)
  dev.off()
  ####
  if(ref=='hg19'|ref=='hg38'){
  print('* Annotate CDBs with CDB conservation. Write outputs:CDB.txt.')
  print('CDB.txt format:#chr start end LRI avgRI conserve_or_not')
  if(ref=='hg19'){
    conserve=read.table(system.file("extdata", "conserved_CDB_hg19.bed", package="RHiCDB"))
  }else{
    conserve=read.table(system.file("extdata", "conserved_CDB_hg19.bed", package="RHiCDB"));
  }
  conserve=cbind(conserve[,1],(conserve[,2]+conserve[,3])/2);
  CDB=cbind(CDB,rep(0,nrow(CDB)));
  for(conchr in 1:23){
  contmp=conserve[conserve[,1]==conchr,2];
  chrpos=which(CDB[,1]==conchr);
  CDBcenter=(CDB[chrpos,2]+CDB[chrpos,3])/2;
  X=as.matrix(contmp)
  Y=as.matrix(CDBcenter)
  dist=outer(seq(length(X)),seq(length(Y)), Vectorize(function(i, j) abs(X[i]-Y[j])))
  pos = apply(dist,2,min) <40000+resolution/2
  CDB[chrpos,6]=pos 
  }
  }else{
  disp('* CDB annotation step is only for hg19 and hg38! Write outputs:CDB.txt.')  
  disp('CDB.txt format:#chr start end LRI avgRI')
  }
  ###
  write.table(CDB,paste(outdir,'/CDB.txt',sep=""),sep="\t",col.names=F,row.names=F,quote=F)
  }
################################################################################################################
  
  ##############################  differential CDB Detection(intersection) ######################################
  
  ##################################################################################################################
}else if(FLAG==2){
    hicfile1=hicfile[[1]][1];
    hicfile2=hicfile[[2]][1];
    # sample 1
    print(paste('A. Processing the first sample ',hicfile1));
    tmp=list.files(path=hicfile1,pattern='chr.*\\.matrix')
    totalchr=length(tmp)
    allpeaks=as.data.frame(matrix(numeric(0),ncol=4))
    allLRI=as.data.frame(matrix(numeric(0),ncol=4))
    allaRI=as.data.frame(matrix(numeric(0),ncol=4))
    print('* Detect all local maximum peaks on each chromosome.')
    
    for( i in 1:totalchr){
      print(paste('Processing chr',i,sep=""))
      hicmap=paste(hicfile1,'/chr',i,'.matrix',sep="")
      out=chrpeaks(hicmap,resolution,chrsizes[i],mind,wd,wdsize,i)
      allpeaks=rbind(allpeaks,out$localmax)
      allLRI=rbind(allLRI,out$LRI)
      allaRI=rbind(allaRI,out$aRI)
    }
    peakidx=order(allpeaks[,4],decreasing = TRUE)
    allpeaks=allpeaks[peakidx,]
    print(' ')
    print('* Write outputs: localmax1.txt, allaRI1.txt, allLRI1.txt')
    write.table(allpeaks,paste(outdir,'/localmax1.txt',sep=""),sep="\t",col.names=F,row.names=F,quote=F)
    write.table(allaRI,paste(outdir,'/allaRI1.txt',sep=""),sep="\t",col.names=F,row.names=F,quote=F)
    write.table(allLRI,paste(outdir,'/allLRI1.txt',sep=""),sep="\t",col.names=F,row.names=F,quote=F)
    print(' ')
    print('* Decide cut off using enrichmentscore. Write outputs:EScurve.png.')
    out=EScutoff(allpeaks,motiffile,outdir)
    CDB1=out$boundary
    score=out$score
    pdf(paste(outdir,'/EScurve1.pdf',sep=""))
    plot(score)
    dev.off()

    
    # sample 2
    print(paste('B. Processing the second sample ',hicfile2));
    tmp=list.files(path=hicfile2,pattern='chr.*\\.matrix')
    totalchr=length(tmp)
    allpeaks=as.data.frame(matrix(numeric(0),ncol=4))
    allLRI=as.data.frame(matrix(numeric(0),ncol=4))
    allaRI=as.data.frame(matrix(numeric(0),ncol=4))
    print('* Detect all local maximum peaks on each chromosome.')
    for( i in 1:totalchr){
      print(paste('Processing chr',i,sep=""))
      hicmap=paste(hicfile2,'/chr',i,'.matrix',sep="")
      out=chrpeaks(hicmap,resolution,chrsizes[i],mind,wd,wdsize,i)
      allpeaks=rbind(allpeaks,out$localmax)
      allLRI=rbind(allLRI,out$LRI)
      allaRI=rbind(allaRI,out$aRI)
    }
    peakidx=order(allpeaks[,4],decreasing = TRUE)
    allpeaks=allpeaks[peakidx,]
    print(' ')
    print('* Write outputs: localmax2.txt, allaRI2.txt, allLRI2.txt')
    write.table(allpeaks,paste(outdir,'/localmax2.txt',sep=""),sep="\t",col.names=F,row.names=F,quote=F)
    write.table(allaRI,paste(outdir,'/allaRI2.txt',sep=""),sep="\t",col.names=F,row.names=F,quote=F)
    write.table(allLRI,paste(outdir,'/allLRI2.txt',sep=""),sep="\t",col.names=F,row.names=F,quote=F)
      print(' ')
      print('* Decide cut off using enrichmentscore. Write outputs:EScurve.png.')
      out=EScutoff(allpeaks,motiffile,outdir)
      CDB2=out$boundary
      score=out$score
      pdf(paste(outdir,'/EScurve2.pdf',sep=""))
      plot(score)
      dev.off()
       ###  compare two samples
      print('')
      print('C.Find differential CDBs.')
      if(ref=='hg19'|ref=='hg38'){
        print('* Annotate differential CDBs with CDB conservation. Write outputs:CDB1.txt,CDB2.txt.')
        print('CDB.txt format:#chr start end LRI avgRI conserve_or_not consistent_or_differential')
        if(ref=='hg19'){
          conserve=read.table(system.file("extdata", "conserved_CDB_hg19.bed", package="RHiCDB"))
        }else{
          conserve=read.table(system.file("extdata", "conserved_CDB_hg38.bed", package="RHiCDB"))
        }
        conserve=cbind(conserve[,1],(conserve[,2]+conserve[,3])/2)
        CDB1=cbind(CDB1,rep(0,nrow(CDB1)),rep(0,nrow(CDB1)))
        CDB2=cbind(CDB2,rep(0,nrow(CDB2)),rep(0,nrow(CDB2)))
        for(conchr in 1:23){
          contmp=conserve[conserve[,1]==conchr,2]
          #conserve  or not  
          chrpos1=which(CDB1[,1]==conchr)
          CDBcenter1=(CDB1[chrpos1,2]+CDB1[chrpos1,3])/2
          pos1=apply(mydist(contmp,CDBcenter1),2,min) < 40000+resolution/2
          CDB1[chrpos1,6]=as.numeric(pos1)
          
          chrpos2=which(CDB2[,1]==conchr)
          CDBcenter2=(CDB2[chrpos2,2]+CDB2[chrpos2,3])/2
          pos2=apply(mydist(contmp,CDBcenter2),2,min)< 40000+resolution/2 
          CDB2[chrpos2,6]=as.numeric(pos2)
          #differential 
          pos=apply(mydist(CDBcenter2[pos2],CDBcenter1[pos1]),2,min)<= 40000  
          CDB1[chrpos1[pos1],7]=as.numeric(pos)
          pos=apply(mydist(CDBcenter1[pos1],CDBcenter2[pos2]),2,min)<= 40000 
          CDB2[chrpos2[pos2],7]=as.numeric(pos)
          
          pos=apply(mydist(CDBcenter2[!pos2],CDBcenter1[!pos1]),2,min)<= resolution+resolution/2  
          CDB1[chrpos1[!pos1],7]=as.numeric(pos)
          pos=apply(mydist(CDBcenter1[!pos1],CDBcenter2[!pos2]),2,min)<= resolution+resolution/2 
          CDB2[chrpos2[!pos2],7]=as.numeric(pos)  
          
        }
      }else{ 
        print('CDB annotation step is only for hg19 and hg38!')
        print('Write outputs without conservation annotation: CDB1.txt, CDB2.txt')      
        print('CDB.txt format:#chr start end LRI avgRI consistent_or_differential')
        CDB1=cbind(CDB1,rep(0,nrow(CDB1)))
        CDB2=cbind(CDB2,rep(0,nrow(CDB2)))
        for(conchr in 1:23){
          #differential 
          chrpos1=which(CDB1[,1]==conchr)
          CDBcenter1=(CDB1[chrpos1,2]+CDB1[chrpos1,3])/2
          chrpos2=which(CDB2[,1]==conchr)
          CDBcenter2=(CDB2[chrpos2,2]+CDB2[chrpos2,3])/2
          
          pos=apply(mydist(CDBcenter2,CDBcenter1),2,min)<= resolution/2+resolution  
          CDB1[chrpos1,6]=as.numeric(pos)
          pos=apply(mydist(CDBcenter1,CDBcenter2),2,min)<= resolution/2+resolution
          CDB2[chrpos2,6]=as.numeric(pos)              
        }
      }
      write.table(CDB1,paste(outdir,'/CDB1.txt',sep=""),sep="\t",row.names=F,col.names=F,quote=F)
      write.table(CDB2,paste(outdir,'/CDB2.txt',sep=""),sep="\t",row.names=F,col.names=F,quote=F)
}else{
  ################################################################################################################
  
  ##############################  differential CDB Detection (test) ######################################
  
  ##################################################################################################################
  n1=length(hicfile[[1]])
  n2=length(hicfile[[2]])
  tmp=list.files(path=hicfile[[1]][1],pattern='chr.*\\.matrix')
  totalchr=length(tmp)
  inni=as.data.frame(matrix(numeric(0),ncol=5))
  allpeaks=list(inni,inni)
  aRI=as.data.frame(matrix(numeric(0),ncol=3+n1+n2))
  size_merge=as.data.frame(matrix(numeric(0),ncol=2))
  print('*A. Calculate the aRI score of each Hi-C replicate.')
  for( i in 1:totalchr){
    print(paste('Processing chr',i,sep=""))
    out=chrpeaks_diff(hicfile,resolution,chrsizes[i],mind,wd,wdsize,i)
    allpeaks[[1]]=rbind(allpeaks[[1]],out$resultA$localmax)
    allpeaks[[2]]=rbind(allpeaks[[2]],out$resultB$localmax)
    size_merge=rbind(size_merge,cbind(out$resultA$aRI[,4],out$resultB$aRI[,4]))
    Ns=ceiling(chrsizes[i]/resolution)-1
    loci=cbind(rep(i,Ns,1),(1:Ns)*resolution-resolution/2,c((1:(Ns-1))*resolution+resolution/2,chrsizes[i]))
    aRI=rbind(aRI,cbind(loci,out$aRI[1:Ns,]))
  }
  
  names(aRI)=c("chr","start","end",NULL)
  print('*B. Get CDBs in the first condition.')
  peakidx=order(allpeaks[[1]][,4],decreasing = TRUE)
  allpeaks[[1]]=allpeaks[[1]][peakidx,]
  print(' ')
  print('* Write outputs: localmax.txt.')
  write.table(allpeaks[[1]],paste(outdir,'/localmax1.txt',sep=""),sep="\t",col.names=F,row.names=F,quote=F)
  print(' ')
  print('* Decide cut off using enrichmentscore. Write outputs:EScurve1.png.')
  out=EScutoff(allpeaks[[1]],motiffile,outdir)
  CDB1=out$boundary
  names(CDB1)=c("chr","start","end","LRI","aRI")
  score=out$score
  pdf(paste(outdir,'/EScurve1.pdf',sep=""))
  plot(score)
  dev.off()
  #
  print('*C. Get CDBs in the second condition.')
  allpeaks[[2]][,4]=allpeaks[[2]][,4]/mean(allpeaks[[2]][,4])*mean(allpeaks[[1]][,4])
  peakidx=order(allpeaks[[2]][,4],decreasing = TRUE)
  allpeaks[[2]]=allpeaks[[2]][peakidx,]
  print(' ')
  print('* Write outputs: localmax.txt.')
  write.table(allpeaks[[2]],paste(outdir,'/localmax2.txt',sep=""),sep="\t",col.names=F,row.names=F,quote=F)
  print(' ')
  print('* Decide cut off using enrichmentscore. Write outputs:EScurve2.png.')
  out=EScutoff(allpeaks[[2]],motiffile,outdir)
  CDB2=out$boundary
  names(CDB2)=c("chr","start","end","LRI","aRI")
  score=out$score
  pdf(paste(outdir,'/EScurve2.pdf',sep=""))
  plot(score)
  dev.off()
  #####################################################
  print('*D. Differential CDBs testing.')
  size=colMeans(aRI[,4:(3+n1+n2)],na.rm =TRUE)
  norm=matrix(rep(size,nrow(aRI)),nrow=nrow(aRI),byrow=TRUE)
  aRI[,4:(3+n1+n2)]=aRI[,4:(3+n1+n2)]/norm*mean(aRI[,4],na.rm =TRUE)
  print('* Write outputs: allaRI.txt. genomic locus + normalized aRI score in each replicates' )
  write.table(aRI,paste(outdir,'/allaRI.txt',sep=""),sep="\t",row.names=F,col.names=F,quote=F)
  CDB=merge(CDB1,CDB2,by=c("chr","start"),all=TRUE)
  CDB=merge(CDB,aRI,by=c("chr","start"),all.x=TRUE,all.y=FALSE)
  CDB$end.x=CDB$end
  CDB=CDB[,-which(names(CDB)%in%c("end.y","end"))]
  CDB=CDB[order(as.numeric(CDB$chr),as.numeric(CDB$start)), ]
  idx=as.data.frame(matrix(numeric(0),ncol=3,nrow(CDB)))
  names(idx)=c("distance","A","B")
  idx$distance=c(FALSE,CDB[2:nrow(CDB),2]==CDB[1:(nrow(CDB)-1),3])
  idx$A=is.na(CDB$aRI.y)
  idx$B=is.na(CDB$aRI.x)
  pos=which(idx$distance)
  CDB_other=CDB[-c(pos,pos-1),]
  CDB_tmp=cbind(CDB[pos-1,4:7],CDB[pos,4:7])
  CDB_tmp[is.na(CDB_tmp[,1]),1:2]=CDB_tmp[is.na(CDB_tmp[,1]),5:6]
  CDB_tmp[is.na(CDB_tmp[,3]),3:4]=CDB_tmp[is.na(CDB_tmp[,3]),7:8]
  CDB_overlap=as.data.frame(matrix(numeric(0),ncol=ncol(CDB),length(pos)))
  names(CDB_overlap)=names(CDB)
  CDB_overlap$chr=CDB[pos-1,]$chr
  CDB_overlap[,2:3]=(CDB[pos-1,2:3]+CDB[pos,2:3])/2
  CDB_overlap[,4:7]=CDB_tmp[,1:4]
  CDB_overlap[,8:(7+n1+n2)]=(CDB[pos-1,8:(7+n1+n2)]+CDB[pos,8:(7+n1+n2)])/2
  CDB=rbind(CDB_other,CDB_overlap)
  
  #size=colMeans(CDB[,8:(7+n1+n2)],na.rm =TRUE)
  #norm=matrix(rep(size,nrow(CDB)),nrow=nrow(CDB),byrow=TRUE)
  #CDB[,8:(7+n1+n2)]=CDB[,8:(7+n1+n2)]/norm*mean(CDB[,8],na.rm =TRUE)
  ##
  test<-function(x,n1,n2){
    if(length(which(is.na(x)))==0){
      ttest=t.test(as.numeric(x[1:n1]),as.numeric(x[(1+n1):(n1+n2)]))
      return(ttest$p.value)}else{
        return(1)
      }
  }
  CDB$diff=rowSums(CDB[,(8+n1):(7+n1+n2)])-rowSums(CDB[,8:(7+n1)])
  CDB$p.value=c(apply(CDB[,8:(7+n1+n2)],1,test,n1=n1,n2=n2))
  print('* Write outputs: CDB_merged.txt. genomic locus + score in merged data + normalized aRI score in each replicates +difference +pvalue')
  write.table(CDB,paste(outdir,'/CDB_merged.txt',sep=""),sep="\t",row.names=F,col.names=F,quote=F)
  
 CDB_A=CDB[which((CDB$p.value<0.05&CDB$diff< -quantile(abs(CDB$diff),0.5,na.rm=TRUE)| CDB$diff< quantile(CDB$diff,0.1,na.rm=TRUE)) & is.na(CDB[,6])),]
 CDB_B=CDB[which((CDB$p.value<0.05&CDB$diff> quantile(abs(CDB$diff),0.5,na.rm=TRUE)| CDB$diff >quantile(CDB$diff,0.9,na.rm=TRUE) )&is.na(CDB[,4])),] 
  
  ###annotation
  if(ref=='hg19'|ref=='hg38'){
    print('* Annotate differential CDBs with CDB conservation. Write outputs:CDB1.txt,CDB2.txt.')
    print('CDB.txt format:#chr start end LRI avgRI conserve_or_not consistent_or_differential')
    if(ref=='hg19'){
      conserve=read.table(system.file("extdata", "conserved_CDB_hg19.bed", package="RHiCDB"))
    }else{
      conserve=read.table(system.file("extdata", "conserved_CDB_hg38.bed", package="RHiCDB"))
    }
    conserve=cbind(conserve[,1],(conserve[,2]+conserve[,3])/2)
    CDB1=cbind(CDB1,rep(0,nrow(CDB1)),rep(0,nrow(CDB1)))
    CDB2=cbind(CDB2,rep(0,nrow(CDB2)),rep(0,nrow(CDB2)))
    for(conchr in 1:23){
      contmp=conserve[conserve[,1]==conchr,2]
      #conserve  or not  
      chrpos1=which(CDB1[,1]==conchr)
      CDBcenter1=(CDB1[chrpos1,2]+CDB1[chrpos1,3])/2
      pos1=apply(mydist(contmp,CDBcenter1),2,min) < 40000+resolution/2
      CDB1[chrpos1,6]=as.numeric(pos1)
      
      chrpos2=which(CDB2[,1]==conchr)
      CDBcenter2=(CDB2[chrpos2,2]+CDB2[chrpos2,3])/2
      pos2=apply(mydist(contmp,CDBcenter2),2,min)< 40000+resolution/2 
      CDB2[chrpos2,6]=as.numeric(pos2)
      #differential 
      chrpos3=which(CDB_A[,1]==conchr)
      CDBcenter3=(CDB_A[chrpos3,2]+CDB_A[chrpos3,3])/2
      chrpos4=which(CDB_B[,1]==conchr)
      CDBcenter4=(CDB_B[chrpos4,2]+CDB_B[chrpos4,3])/2
      
      pos=apply(mydist(CDBcenter3,CDBcenter1),2,min)<= resolution/2+1 
      CDB1[chrpos1,7]=1-as.numeric(pos)
      pos=apply(mydist(CDBcenter4,CDBcenter2),2,min)<= resolution/2+1 
      CDB2[chrpos2,7]=1-as.numeric(pos)        
    }
  }else{ 
    print('CDB annotation step is only for hg19 and hg38!')
    print('Write outputs without conservation annotation: CDB1.txt, CDB2.txt')      
    print('CDB.txt format:#chr start end LRI avgRI consistent_or_differential')
    CDB1=cbind(CDB1,rep(0,nrow(CDB1)))
    CDB2=cbind(CDB2,rep(0,nrow(CDB2)))
    for(conchr in 1:23){
      #differential 
      chrpos1=which(CDB1[,1]==conchr)
      CDBcenter1=(CDB1[chrpos1,2]+CDB1[chrpos1,3])/2
      chrpos2=which(CDB2[,1]==conchr)
      CDBcenter2=(CDB2[chrpos2,2]+CDB2[chrpos2,3])/2
      chrpos3=which(CDB_A[,1]==conchr)
      CDBcenter3=(CDB_A[chrpos3,2]+CDB_A[chrpos3,3])/2
      chrpos4=which(CDB_B[,1]==conchr)
      CDBcenter4=(CDB_B[chrpos4,2]+CDB_B[chrpos4,3])/2
      
      pos=apply(mydist(CDBcenter3,CDBcenter1),2,min)<= resolution/2+1  
      CDB1[chrpos1,6]=1-as.numeric(pos)
      pos=apply(mydist(CDBcenter4,CDBcenter2),2,min)<= resolution/2+1 
      CDB2[chrpos2,6]=1-as.numeric(pos)              
    }
  }
  write.table(CDB1,paste(outdir,'/CDB1.txt',sep=""),sep="\t",row.names=F,col.names=F,quote=F)
  write.table(CDB2,paste(outdir,'/CDB2.txt',sep=""),sep="\t",row.names=F,col.names=F,quote=F)
  
  
}
  
  }

