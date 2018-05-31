# RHiCDB – user guide
## Overview
&emsp;&emsp;RHiCDB is an open-source R package based on HiCDB methods that detects the contact domain boundaries (CDBs) from Hi-C contact matrix. RHiCDB function takes raw or normalized contact matrix and outputs conservation annotated CDBs or differential CDBs. visHiCDB function takes raw or normalized contact matrix and HiCDB results and outputs visualization of CDBs on single Hi-C map or differential CDBs on two Hi-C maps. HiCDB are also implemented as [MATLAB version](https://github.com/ChenFengling/HiCDB).  
&emsp;&emsp;Here is the general features of HiCDB.  
<div align=center><img width="420" height="300" src="https://github.com/ChenFengling/RHiCDB/blob/master/images/pipeline.png"/></div>

&emsp;&emsp;Here is the general steps of how we detect CDBs.  
<div align=center><img width="580" height="380" src="https://github.com/ChenFengling/RHiCDB/blob/master/images/pipe.png"/></div>

### Requirements and install
Download RHiCDB_1.0.tar.gz and install in R. RHiCDB depends on pracma,limma,Matrix,gridExtra,rasterVis and lattice.
```R
install.packages("RHiCDB_1.0.tar.gz")
```
### Quick start
Unzip the testdata.tar.gz, you will  find the dense format Hi-C data of hESC (Doxin et al.). 
```shell
tar -zxvf testdata.tar.gz
```
##### Run RHiCDB to get the CDBs.
```R
library('RHiCDB')
hicfile='/home/fchen/h1_rep1/'
resolution=40000
chrsizes='hg19'
outdir='/home/fchen/h1_rep1/'
RHiCDB(hicfile,resolution,chrsizes,ref='hg19',outdir=outdir)
```
This will take the  intra-chromosome matrix ('chr1.matrix',...,'chr23.matrix') in '/home/fchen/h1_rep1/' as input and set the resolution as 40000,chrsizes as 'hg19', the CTCF motif ref as 'hg19' and output the contact domain boundaries.   
##### Run visHiCDB to display the region chr17:67100000-71100000.
```R
hicfile='/home/fchen/h1_rep1/chr17.matrix'
resolution=40000
outdir='/home/fchen/h1_rep1'
CDBfile='/home/fchen/h1_rep1/CDB.txt'
chr=17
startloc=67100000
endloc=71100000
visHiCDB(hicfile,CDBfile,resolution,chr,startloc,endloc,outdir)
```
You will get this output in 17_67100000_71100000_HiCmap.pdf. The dot is CDB detected(dark blue:consistently detected CDBs; light blue:other CDBs)  
<div align=center><img width="500" height="500" src="https://github.com/ChenFengling/RHiCDB/blob/master/images/17_67100000_71100000_HiCmap.png"/></div>

#####  Get .bed file
As "chrX" is named as "chr23" and as "23" in the output CDB.txt file. You could use the following shell code to change CDB.txt into .bed file.    
```shell
awk -v OFS="\t" '{ print "chr"$1,$2,$3,$4,$5}' CDB.txt >CDB.bed
sed -i  's/chr23/chrX/g' CDB.bed
```  
## 1. Run RHiCDB 
### Input
**hicfile:** The directory of all intra-chromosome matrix of a sample. The intra-chromosome matrix must be named as "chr+number.matrix" according to the chromosome order like 'chr1.matrix','chr2.matrix',...,'chr23.matrix'. As HiCDB matches "chr\*.matrix" to recognize the Hi-C matrix, avoid to use the "chr\*.matrix" as the name of other files. The intra-chromosome matrix could be in a dense (a NxN matrix) or sparse (a Kx3 table,Rao et al.) format. hicfile should be set as 'SAMPLE_DIR' when option is "singlemap", list('SAMPLE_DIR1','SAMPLE_DIR2') or  list(c(’SAMPLE1_rep1’,’SAMPLE1_rep2’),c(’SAMPLE2_rep1’,’SAMPLE2_rep2’)) when option is ‘comparemap’. This is required.  
Dense format contains	the contact	frequencies	of the Hi-C NxN	matrix.  
Sparse format (Rao et al.) has three fields: i, j, and M_i,j. (i and j are written as the left edge of the bin at a given resolution; for example, at 10 kb resolution, the entry corresponding to the first row and tenth column of the matrix would correspond to M_i,j, where i=0, j=90000). As the Hi-C matrix is symmetric, only the upper triangle of the matrix is saved in sparse format. An example is as following:     

| | | |
|-|-|-|
|  50000   |  50000    |  1.0 |
|  60000   |  60000    |  1.0 |
|  540000  |  560000   |  1.0 |
|  560000  | 560000    | 59.0 |
|  560000  |  570000   |  1.0 |
|  560000  |  600000   |  1.0 |
|  560000  |  700000   |  1.0 |
|  690000  |  710000   |  1.0 |
|  700000  |  710000   |  1.0 |
|  710000  | 710000    | 66.0 |

**resolution:** resolution of Hi-C matrix. This is required.  
**chrsizes:** Ordered chromosome sizes of the genome. Optional setting is ‘hg19’, ‘hg38’, ‘mm9’, ‘mm10’ or any other chromosome size files which can be generated following the instructions in annotation/README.md. This is required.  
**ref:** ref should be set when you want to get a cutoff using a CTCF motif or the option is 'comparemap'. Optional ref is ‘hg19’, ‘hg38’, ‘mm9’, ‘mm10’ or any other custom motif locus files which can be generated from instructions in annotation/README.md. Only ‘hg19’ and ‘hg38’ can be annotated with conservation.  
### Examples
#### 1. Output all the local maximum peaks and let customers to decide the cutoff.  
```R
RHiCDB('sample1/',10000,chrsizes='custom_chrsizes.txt');
RHiCDB('sample1/',10000,chrsizes='custom_chrsizes.txt',outdir='sample1/outputs/');
```
#### 2. Use GSEA-like methods to decide the cutoff . 
```R
RHiCDB('sample1/',10000,chrsizes='hg19',ref='hg19');
RHiCDB('sample1/',10000,chrsizes='custom_chrsizes.txt',ref='custom_motiflocs.txt')
```
#### 3. To detect differential CDBs
```R
RHiCDB(list('sample1','sample2'),10000,'hg19',ref='hg19');
RHiCDB(list(c("sample1_rep1","sample1_rep2"),c("sample2_rep1","sample2_rep2")),10000,'hg19',ref='hg19');
```
### Output(s)
**1.CDB.txt:**

| chr | start|  end|    LRI|    avgRI   | conserve_or_not |     consistent_or_differential |
|:------:|:----:|:------:|:----:|:------:|:----:|:----:|
| 19    | 53100000      | 53140000      | 0.394707211   | 0.647392804 | 0 |     1 |
| 16    | 5060000         | 5100000       | 0.342727704 | 0.663101081 | 1 |     1 |
| 19    | 19620000  | 19660000  |       0.329837698 | 0.609237673 |     1 |     0 |

**2. localmax.txt:** all the local maximum peaks detected before cutoff decision. User can decide custom CDB cutoff upon this file.  
**3. EScurve.png:** CTCF motif enrichment on ranked local maximum peaks.  
These output files can be found in custom output directory or default directory namely the directory of the first sample.      
**4. aRI.txt:** average RI score for each  genomic bin.  
**5. LRI.txt:** LRI score for each genomic bin.
## 2. Run visHiCDB
### Input
**hicfile:** Hi-C matrix of the intersested chromosome.
**CDBfile:** CDBfile sould be a cell array storing the CDB location. The CDB files should be formatted as the output files of function HiCDB.  
**resolution:** resolution of Hi-C map.  
**chr,startloc,endloc:** observation locus on Hi-C map.  
### Examples
#### 1.Show CDB on single Hi-C map
```R
visHiCDB('sample1/chr18.matrix','CDB1.txt',40000,18,25000000,31150000)
```
#### 2. Show differential CDBs on Hi-C maps 
```R
visHiCDB(list('sample1/chr18.matrix','sample2/chr18.mat'),list('CDB1.txt','CDB2.txt'),40000,18,25000000,31150000)
```
### Output(s)
**HiCmap.pdf：** a pdf containing figure showing CDBs on single Hi-C map or different kinds of CDBs between two Hi-C maps.  
These output files can be found in custom output directory or default directory namely the directory of the first sample. 


