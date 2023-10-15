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

#install.packages("readxl")
library("readxl")
##################################################################################################

args<-commandArgs(TRUE)
window_list = list(as.numeric(50000))   ##in bp
species1<-'hg38'			## UCSC stle abbreviation (hg19, gorGor3, etc)	
species2<-'gorGor6'		
## either "all" or UCSC style notation (chr1, chr2, chrX, etc)
chrom_list=list('chrY')
chrom_list2='all'
#percent_list = list(as.numeric(100), as.numeric(98), as.numeric(96), as.numeric(96.6), as.numeric(90))  #it means the CNE identification threshold
percent_list = list(as.numeric(100))


for(window_size in window_list)
{
      for(percent in percent_list)
      {
            for(chrom in chrom_list)
            {   
                #"/bucket/.mabuya/MillerU/Abrar/PCS/chrY_MSY_coordinates.xlsx"
                identical_seqs.gr<-makeGRangesFromDataFrame(read_excel("/bucket/.mabuya/MillerU/Abrar/PCS/chrY_MSY_coordinates.xlsx"))
                 
                #print(identical_seqs.gr)
                
                ####tile the genome
                
                genome_info<-read.table(paste0("/bucket/MillerU/Abrar/bigZips/", species1, ".chrom.sizes"))
                
                genome_seqlengths<-genome_info$V2
                names(genome_seqlengths)<-genome_info$V1
                
                tiled_genome<-tileGenome(genome_seqlengths, tilewidth=window_size, cut.last.tile.in.chrom=T)
                
                if(chrom=="all"){
                	tiled_chr<-tiled_genome
                } else {
                	tiled_chr<-tiled_genome[seqnames(tiled_genome)==chrom]
                }
                
                ####get distribution in each window and calculate kurtosis
                
                tileHits<-findOverlaps(tiled_chr, identical_seqs.gr)
                tileSplit<-split(subjectHits(tileHits), queryHits(tileHits))
                #print(identical_seqs.gr)
                print(tiled_chr)
                write.table( x = data.frame(tileHits), file = paste0("/bucket/.mabuya/MillerU/Abrar/PCS/hg38-gorGor6/MSY_intersected_",chrom,".csv"), sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE )
            }
      }
}



#gc()