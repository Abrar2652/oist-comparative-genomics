#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

#BiocManager::install("GenomicRanges")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

#BiocManager::install("rtracklayer")
#install.packages('rtracklayer')
#install.packages('data.table')
#install.packages("remotes")
#library(remotes)
#install.github('kloke/npsm')
#install.packages('devtools')
#devtools::install_github("CshlSiepelLab/RPHAST")
#install.packages('rphast',"/home/m/md-jahin/R/x86_64-redhat-linux-gnu-library/3.6",'http://ftp.ussg.iu.edu/CRAN')
library(GenomicRanges)
library(rtracklayer)
library(data.table)
library(npsm)
library(rphast)
library(e1071)  

##################################################################################################

args<-commandArgs(TRUE)

window_size<-as.numeric(30000) 	##in bp

species1<-'hg19'			## UCSC stle abbreviation (hg38, calJac4, etc)	
#species2<-list('calJac4'=30,'danRer10'=8,'equCab3'=12,'felCat9'=10,'galGal6'=8,'ornAna2'=7,'oryCun2'=10,'oviAri4'=11,'petMar3'=9,'susScr11'=10,'tarSyr2'=14,'thaSir1'=8,'triMan1'=10,'xenTro10'=10)
## either "all" or UCSC style notation (chr1, chr2, chrX, etc)
species2<-list('mm10'=1)
#chrom_list=list('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY')
chrom_list=list('all')
percent=as.numeric(100) #it means the CNE identification threshold

for(i in seq_along(species2))
{
    outDir<-paste0("/bucket/MillerU/Abrar/hg19-", names(species2)[i], "_outputs/hg19-", names(species2)[i], "_filtered_pcs/")
    
    if(!dir.exists(file.path(outDir))){
    	dir.create(file.path(outDir))
    }
    
    for(chrom in chrom_list)
    {    
        ## read in the identical sequence locations (tsv files) - output of get_identical_seq_locations.R
        # get_identical_seq_loc_hg19-calJac4_L20    
        # filtered_pcs_hg19-calJac4_filtered_L10    
        identical_seqs_fn<-paste0("/bucket/MillerU/Abrar/hg19-", names(species2)[i], "_outputs/get_identical_seq_loc_hg19-", names(species2)[i], "/", chrom)
        print(identical_seqs_fn)
        if(chrom == "all"){
          identical_seqs_files<-dir(identical_seqs_fn, pattern=paste0(species1, "_", names(species2)[i], ".*all*"), full.names=T)
        } else {
          identical_seqs_files<-dir(identical_seqs_fn, pattern=paste0(species1, "_", names(species2)[i], ".*chr*"), full.names=T)
        }
    
        identical_seqs.dt<-rbindlist(lapply(identical_seqs_files, fread))
        	
        identical_seqs.gr<-GRanges(seqnames=identical_seqs.dt$first.seqnames, 
        		   ranges=IRanges(start=identical_seqs.dt$first.start,
                end=identical_seqs.dt$first.end))
        
        if(chrom != "all"){
        	identical_seqs.gr<-subset(identical_seqs.gr, seqnames(identical_seqs.gr)==chrom)
        }
    
        print(identical_seqs.gr)
        ### read in the filter regions - usually a list of exon and repeat coordinates
        filter_regions <- import.bed(dir(paste0("/bucket/MillerU/Abrar/HCNE_hg19_mm10/100/"), pattern=".*.bed", full.names=T)) 

        ##Filter regions for reference genome - repeats and exons
        filtered_identical_seqs.gr<-setdiff(identical_seqs.gr, filter_regions)
        
        write.table( x = data.frame(filtered_identical_seqs.gr), file = paste0(outDir, species1, "_", species2, "_", chrom, "_", i, "_identical_seqs.csv"), sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE )
        #print("tiled pcs calculated")
    }
}

#gc()