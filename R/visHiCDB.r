visHiCDB=function(hicfile,CDBfile,resolution,chr,startloc,endloc,outdir){
  #'@title   Visualization of CDBs or differential CDBs on Hi-C maps
  #' 
  #'@description  visHiCDB uses Hi-C raw contact matrix and CDBs as input and outputs figures
  #'   of CDBs or differential CDBs on Hi-C maps
  #' @details    This function outputs a pdf figure showing CDBs on a Hi-C
  #'       map on desired locus. Conserved CDBs are marked as dark blue.
  #'@param hicfile hicfile is the file names of  the intra-chromosome matrixes.
  #'      The intra-chromosome matrix could be in a dense (a NxN matrix) or sparse (a Kx3 table,Rao et al.) format. 
  #'   Show CDBs on one sample,set hicfile as 'SAMPLE_File'.  
  #'   Show  differebtial CDBs in two samples, set hicfile as list('SAMPLE1_FILE','SAMPLE2_FILE').
  #'@param   CDBfile  CDBfile should be is file name of the CDB files.
  #'         Show CDBs on one sample,set CDBfile as 'SAMPLE_CDB'.  
  #'       Show  differebtial CDBs in two samples, set CDBfile as list('SAMPLE1_CDB','SAMPLE2_CDB').
  #'    The CDB file should be formated as the output file of HiCDB.
  #' @param resolution  resolution of Hi-C matrix. This is required.
  #'@param  chr,startloc,endloc   numeric observation locus on Hi-C map.This is required.
  #' @param  outdir   The output direction. Default will be the directory of  the first sample.       
  #'@author
  #'   Implemented by Fengling Chen
  #' 
  #'   Any suggestions and remarks might be addressed to Fengling Chen:cfl15@mails.tsinghua.edu.cn 
  #' @export  
  #' @examples
  #'      ## Show CDB on single Hi-C map
  #'      visHiCDB('sample1/chr18.matrix','CDB1.txt',40000,18,25000000,31150000)
  #'      ## Show differential CDBs
  #'      visHiCDB(list('sample1/chr18.matrix','sample2/chr18.matrix'),
  #'         +list('CDB1.txt','CDB2.txt'),40000,18,25000000,31150000)

  .loadPackage("pracma") 
  .loadPackage("lattice")
  .loadPackage("rasterVis")
  .loadPackage("gridExtra")
  
  ################################################################################################################
  
  ##############################  check input ######################################
  
  ##################################################################################################################
  if (length(hicfile)==1){
    print("One sample is provided, show its heatmap!")
          FLAG=1
    if (missing(outdir)){
      outdir=dirname(hicfile[1])
    }else if(!dir.exists(outdir)){
       dir.create(outdir,recursive=TRUE)
    }  
  }else if(length(hicfile==2)){
          FLAG=2
      print('Two sample is provided,show comparing heatmap!')
    if (missing(outdir)){
      outdir=dirname(hicfile[[1]][1])
    }else if(!dir.exists(outdir)){
       dir.create(outdir,recursive=TRUE)
    }

  }else{
    stop("You should provided samples following the documents!")
  }
   print(paste('Outdir:',outdir))
  
  if(length(hicfile)!=length(CDBfile)){
    stop('unequal hicfile and CDB file length.')
  }
  
  if (is.null(resolution)){
    stop("resolution is required!")
  }
  

  
  tX=round(startloc/resolution):round(endloc/resolution)
  ################################################################################################################
  
  ##############################  CDB detection ######################################
  
  ##################################################################################################################
  if(FLAG==1){
    im=read2dense(hicfile,resolution=resolution)
    out=KRnorm(im)
    imnew=out$imKR
    
    CDB=read.table(CDBfile,stringsAsFactors=FALSE)
    CDB=CDB[CDB[,1]==chr&CDB[,2]>=(tX[1]-1)*resolution & CDB[,3]<=tX[length(tX)]*resolution,] 
    bound=(CDB[,2]+CDB[,3])/resolution/2-tX[1]+1+0.5
    obj <- data.frame(x=bound, y=bound)
    coordinates(obj) <- ~x+y
    mat=log(1+imnew[tX,tX])
    colnames(mat)=NULL
    rownames(mat)=NULL
    #p=levelplot(mat,par.settings = YlOrRdTheme,margin=FALSE,xlab = list(label = ""),ylab = list(label = ""))
    pdf(paste(outdir,"/",chr,'_',startloc,'_',endloc,'_HiCmap.pdf',sep=""),7,7)
    if(ncol(CDB)>5){
      #p=p+layer(sp.points(obj,col=ifelse(CDB[,6]==1,'#002060','#5B9BD5'),pch=19,cex=0.8)) 
     p=levelplot(mat,par.settings = YlOrRdTheme,margin=FALSE,xlab = list(label = ""),ylab = list(label = ""),panel=function(...){
        panel.levelplot(...)
        sp.points(obj,col=ifelse(CDB[,6]==1,'#002060','#5B9BD5'),pch=19,cex=0.8)
      })
      }else{
        p=levelplot(mat,par.settings = YlOrRdTheme,margin=FALSE,xlab = list(label = ""),ylab = list(label = ""),panel=function(...){
          panel.levelplot(...)
          sp.points(obj,col='#000000',pch=19,cex=0.8)
        })
    }
    print(p)
    dev.off()
    

    ################################################################################################################
    
    ##############################  differential CDB Detection ######################################
    
    ##################################################################################################################
  }else{
    pdf(paste(outdir,"/",chr,'_',startloc,'_',endloc,'_HiCmap_compare.pdf',sep=""),9,3.9)
          print('*Processing sample 1*')
          im=read2dense(hicfile[[1]][1],resolution=resolution)
          out=KRnorm(im)
          imnew=out$imKR
           depth=sum(imnew)         
          CDB=read.table(CDBfile[[1]][1],stringsAsFactors=FALSE)
          CDB=CDB[CDB[,1]==chr&CDB[,2]>=(tX[1]-1)*resolution & CDB[,3]<=tX[length(tX)]*resolution,] 
          bound=(CDB[,2]+CDB[,3])/resolution/2-tX[1]+1+0.5
          obj1 <- data.frame(x=bound, y=bound)
          coordinates(obj1) <- ~x+y
          mat=log(1+imnew[tX,tX])
          colnames(mat)=NULL
          rownames(mat)=NULL

          if(ncol(CDB)==5){
            p1=levelplot(mat,par.settings = YlOrRdTheme,margin=FALSE,xlab = list(label = ""),ylab = list(label = ""),panel=function(...){
              panel.levelplot(...)
              sp.points(obj1,col='#000000',pch=19,cex=0.5)
            })
          }else if(ncol(CDB)==6){
            p1=levelplot(mat,par.settings = YlOrRdTheme,margin=FALSE,xlab = list(label = ""),ylab = list(label = ""),panel=function(...){
              panel.levelplot(...)
              sp.points(obj1,col=ifelse(CDB[,6]==1,'#002060','#5B9BD5'),pch=19,cex=0.5)
            })
          }else if(ncol(CDB)==7){
             color1=c()
             color1[CDB[,6]&CDB[,7]]='#002060'
             color1[!CDB[,6]&CDB[,7]]='#5B9BD5'
             color1[CDB[,6]&!CDB[,7]]='#ED7D31'
             color1[!CDB[,6]&!CDB[,7]]='#FFFF68'
             p1=levelplot(mat,par.settings = YlOrRdTheme,margin=FALSE,xlab = list(label = ""),ylab = list(label = ""),panel=function(...){
               panel.levelplot(...)
               sp.points(obj1,col=color1,pch=19,cex=0.5)
             })
            }
          

        print('*Processing sample 2*')
        im=read2dense(hicfile[[2]][1],resolution=resolution)
        out=KRnorm(im)
        imnew=out$imKR
        imnew=imnew/sum(imnew)*depth
        CDB=read.table(CDBfile[[2]][1],stringsAsFactors=FALSE)
        CDB=CDB[CDB[,1]==chr&CDB[,2]>=(tX[1]-1)*resolution & CDB[,3]<=tX[length(tX)]*resolution,] 
        bound=(CDB[,2]+CDB[,3])/resolution/2-tX[1]+1+0.5
        obj2 <- data.frame(x=bound, y=bound)
        coordinates(obj2) <- ~x+y
        mat=log(1+imnew[tX,tX])
        colnames(mat)=NULL
        rownames(mat)=NULL
        p2=levelplot(mat,par.settings = YlOrRdTheme,margin=FALSE,xlab = list(label = ""),ylab = list(label = ""))
        
        if(ncol(CDB)==5){
          p2=levelplot(mat,par.settings = YlOrRdTheme,margin=FALSE,xlab = list(label = ""),ylab = list(label = ""),panel=function(...){
            panel.levelplot(...)
            sp.points(obj2,col='#000000',pch=19,cex=0.8)
          })
        }else if(ncol(CDB)==6){
          p2=levelplot(mat,par.settings = YlOrRdTheme,margin=FALSE,xlab = list(label = ""),ylab = list(label = ""),panel=function(...){
            panel.levelplot(...)
            sp.points(obj2,col=ifelse(CDB[,6]==1,'#002060','#5B9BD5'),pch=19,cex=0.5)
          })
        }else if(ncol(CDB)==7){
          color2=c()
          color2[CDB[,6]&CDB[,7]]='#002060'
          color2[!CDB[,6]&CDB[,7]]='#5B9BD5'
          color2[CDB[,6]&!CDB[,7]]='#7030A0'
          color2[!CDB[,6]&!CDB[,7]]='#99FF99'
          p2=levelplot(mat,par.settings = YlOrRdTheme,margin=FALSE,xlab = list(label = ""),ylab = list(label = ""),panel=function(...){
            panel.levelplot(...)
            sp.points(obj2,col=color2,pch=19,cex=0.5)
          })
        }
        grid.arrange(p1=p1,p2=p2,ncol=2)
        dev.off()
        return(list(color1,color2))
}

}
