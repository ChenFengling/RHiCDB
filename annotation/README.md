## How to generate CTCF motif sites file ?  
&emsp;&emsp;These files can be generated using homer.  
&emsp;&emsp;Install homer and scan CTCF motif on certain genome.  
```
##scan motif
perl HOMER_PATH/.//configureHomer.pl -install hg19
scanMotifGenomeWide.pl  HOMER_PATH/motifs/ctcf.motif  hg19 -txt > CTCF_hg19_site.txt

##modify CTCF sites file to perfered format  
##human
sed 's/chr//g' CTCF_hg19_site.txt >tmp
awk '{printf("%s\t%d\n",$1,($2+$3)/2)}' tmp > CTCF_hg19_site.txt
sed -i '/_/d' CTCF_hg19_site.txt
sed -i 's/X/23/g' CTCF_hg19_site.txt
sed -i 's/Y/24/g' CTCF_hg19_site.txt

##mouse
sed 's/chr//g' ctcf.sites.mm10.txt >tmp
awk '{printf("%s\t%d\n",$1,($2+$3)/2)}' tmp > CTCF_mm10_site.txt
sed -i '/_/d' CTCF_mm10_site.txt
sed -i 's/X/20/g' CTCF_mm10_site.txt
sed -i 's/Y/21/g' CTCF_mm10_site.txt
```

## How to generate chromosome size files ?
&emsp;&emsp;The chromosome size files can be downloaded on the UCSC genome browser and be reformatted as the lengths of sorted chromosomes.  
```
cat hg38.chrom.sizes|sed  '/_/d'|sed '/chrM/d'|sort -k1,1 -V -s|cut -f2 >hg38.chrom.sizes.txt
``` 

## How are the conserved CDB files generated ?
&emsp;&emsp;conserved_CDB_hg19.bed was the CDBs detected in at least 16 samples of the 21 Hi-C samples from Schmitt et al(1).

**reference**  
(1)Schmitt, A.D., Hu, M., Jung, I., Xu, Z., Qiu, Y., Tan, C.L., Li, Y., Lin, S., Lin, Y. and Barr, C.L. et al.. (2016) A Compendium of Chromatin Contact Maps Reveals Spatially Active Regions in the Human Genome. CELL REP, 17, 2042-2059.
