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

window_size<-as.numeric(50000) 	##in bp  [10000000, 5000000, 2500000, 1200000, 500000, 100000, 50000, 30000, 15000]

species1<-'hg38'			## UCSC stle abbreviation (hg38, calJac4, etc)	
species2<-list('gorGor6'=1)

#('calJac4'=30,'danRer10'=8,'equCab3'=12,'felCat9'=10,'galGal6'=8,'ornAna2'=7,'oryCun2'=10,'oviAri4'=11,'petMar3'=9,'susScr11'=10,'tarSyr2'=14,'thaSir1'=8,'triMan1'=10,'xenTro10'=10)
## either "all" or UCSC style notation (chr1, chr2, chrX, etc)

chrom_list=list('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY')
percent=as.numeric(100) #it means the CNE identification threshold

for(i in seq_along(species2))
{
    outDir<-paste0("/bucket/MillerU/Abrar/hg38-", names(species2)[i], "_outputs/hg38-", names(species2)[i], "_sliding_pcs_", window_size/1000, "kb/")
    
    if(!dir.exists(file.path(outDir))){
    	dir.create(file.path(outDir))
    }
    
    for(chrom in chrom_list)
    {    
        ## read in the identical sequence locations (tsv files) - output of get_identical_seq_locations.R
        # get_identical_seq_loc_hg38-calJac4_L20    
        # filtered_pcs_hg38-calJac4_filtered_L10    
        #identical_seqs_fn<-paste0("/bucket/MillerU/Abrar/hg38-", names(species2)[i], "_outputs/get_identical_seq_loc_hg38-", names(species2)[i], "_vanilla_L", species2[[i]], "/", chrom)
        identical_seqs_fn<-paste0("/bucket/MillerU/Abrar/hg38-", names(species2)[i], "_outputs/get_identical_seq_loc_hg38-", names(species2)[i], "/", chrom)
        print(identical_seqs_fn)
        if(chrom == "all"){
          identical_seqs_files<-dir(identical_seqs_fn, pattern=paste0(species1, "_", names(species2)[i], ".*all*"), full.names=T)
        } else {
          identical_seqs_files<-dir(identical_seqs_fn, pattern=paste0(species1, "_", names(species2)[i], ".*chr*"), full.names=T)
        }
    
        identical_seqs.dt<-rbindlist(lapply(identical_seqs_files, fread))
        	
        identical_seqs.gr<-GRanges(seqnames=identical_seqs.dt$first.seqnames, 
        		   ranges=IRanges(start=identical_seqs.dt$first.start,
                end=identical_seqs.dt$first.end),
        		   strand="*")
        
        if(chrom != "all"){
        	identical_seqs.gr<-subset(identical_seqs.gr, seqnames(identical_seqs.gr)==chrom)
        }
    
    
        ##Filter regions for reference genome - repeats and exons
        filtered_identical_seqs.gr <- identical_seqs.gr
        
        ####tile the genome
        ##downloaded the .sizes file from https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/
        
        genome_info<-read.table(paste0("/bucket/MillerU/Abrar/bigZips/", species1, ".chrom.sizes"))
        
        gr <- GRanges(
                      seqnames=genome_info$V1,
                      ranges=IRanges(start=1,
                      end=genome_info$V2),
                      strand="*"
                      )
        
        tiled_genome <- unlist(as((slidingWindows(gr, width=50000L, step=1000L)), "GRangesList"))
        print(tiled_genome)
        
        if(chrom=="all"){
        	tiled_chr<-tiled_genome
        } else {
        	tiled_chr<-tiled_genome[seqnames(tiled_genome)==chrom]
        }   
        ####get distribution in each window and calculate kurtosis 
        tileHits<-findOverlaps(tiled_chr, filtered_identical_seqs.gr)
        tileSplit<-split(subjectHits(tileHits), queryHits(tileHits))
        
        pcs_by_tile<-lapply(as.list(1:length(tiled_chr)) , function(x){
                    obj<-tileSplit[[as.character(x)]]
                    if(is.null(obj)){
                      return(0)
                    } else {
                      runs<-filtered_identical_seqs.gr[obj]
                      dist<-width(runs)
                      score<-dist
                      if(is.na(score)) score = 0
                      if(score=="NaN") score = 0
                      return(score)
                    }
                    })
        				
        tiled_chr$score<-pcs_by_tile
        #drop windows that aren't the correct width - the ends of chromosomes usually
        tiled_chr<-tiled_chr[width(tiled_chr)==window_size]
    
        dt_text <- unlist(lapply(tiled_chr$score, paste, collapse=","))   
        writeLines(dt_text, paste0(outDir, chrom, "_", window_size/1000, "kb_min_window_sliding_pcs.txt"))
        print("sliding pcs calculated")
        
        
    }
}

#gc()