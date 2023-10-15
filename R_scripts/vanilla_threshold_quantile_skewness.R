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


##################################################################################################

args<-commandArgs(TRUE)
window_list = list(as.numeric(30000), as.numeric(15000))   #list(as.numeric(30000), as.numeric(15000), as.numeric(120000), as.numeric(60000))   ##in bp
species1<-'hg38'			## UCSC stle abbreviation (hg38, micMur3, etc)	
species2<-list('calJac4'=30,'danRer10'=8,'equCab3'=12,'felCat9'=10,'galGal6'=8,'ornAna2'=7,'oryCun2'=10,'oviAri4'=11,'petMar3'=9,'susScr11'=10,'tarSyr2'=14,'thaSir1'=8,'triMan1'=10,'xenTro10'=10)			
## either "all" or UCSC style notation (chr1, chr2, chrX, etc)
chrom_list=list('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY')
chrom_list2='all'
#percent_list = list(as.numeric(100), as.numeric(98), as.numeric(96), as.numeric(96.6), as.numeric(90))  #it means the CNE identification threshold
percent_list = list(as.numeric(100))

for(i in seq_along(species2))
{
    for(chrom in chrom_list2)
    {    
        ## read in the identical sequence locations (tsv files) - output of get_identical_seq_locations.R
        
        ## get_identical_seq_loc_hg38-micMur3   
        ## filtered_pcs_hg38-micMur3   
        identical_seqs_fn<-paste0("/bucket/MillerU/Abrar/hg38-", names(species2)[i], "_outputs/get_identical_seq_loc_hg38-", names(species2)[i], "_vanilla_L", species2[[i]], "/", chrom)
        
        if(chrom == "all"){
          identical_seqs_files<-dir(identical_seqs_fn, pattern=paste0(species1, "_", names(species2)[i], ".*all*"), full.names=T)
        } else {
          identical_seqs_files<-dir(identical_seqs_fn, pattern=paste0(species1, "_", names(species2)[i], ".*chr*"), full.names=T)
        }
        print(identical_seqs_files)
        identical_seqs.dt<-rbindlist(lapply(identical_seqs_files, fread))
        
        ###get the range of the middle 50% of the distribution of all runs of conserved seqs across the genome - used later in kurtosis calculation    
        denominator_quants<-quantile(identical_seqs.dt$first.width, c(0.75, 0.25))
        denominator<-denominator_quants[1]-denominator_quants[2]
        print(denominator)	
    }
    
    for(window_size in window_list)
    {
          for(percent in percent_list)
          {
                for(chrom in chrom_list)
                {   
                    #outDir<-paste0("/bucket/MillerU/Abrar/hg38-micMur3_outputs/", window_size/1000, "kb/binned_quantile_kurtosis_genome_ave", "_", percent, "pc_500col/") #Gorilla
                    outDir<-paste0("/bucket/MillerU/Abrar/hg38-", names(species2)[i], "_outputs/", window_size/1000, "kb/binned_quantile_skewness_genome_ave", "_", percent, "pc/") #Others
                    
                    if(!dir.exists(file.path(outDir))){
                    	dir.create(file.path(outDir))
                    } 
                    ## read in the identical sequence locations (tsv files) - output of filter_exons_repeats.py   
                    ## get_identical_seq_loc_hg38-micMur3   
                    # filtered_pcs_hg38-micMur3
                    identical_seqs_fn<-paste0("/bucket/MillerU/Abrar/hg38-", names(species2)[i], "_outputs/get_identical_seq_loc_hg38-", names(species2)[i], "_vanilla_L", species2[[i]], "/", chrom)
                    
                    if(chrom == "all"){
                      identical_seqs_files<-dir(identical_seqs_fn, pattern=paste0(species1, "_", names(species2)[i], ".*all*"), full.names=T)
                    } else {
                      identical_seqs_files<-dir(identical_seqs_fn, pattern=paste0(species1, "_", names(species2)[i], ".*chr*"), full.names=T)
                    }
                    print(identical_seqs_files)
                    identical_seqs.dt<-rbindlist(lapply(identical_seqs_files, fread))
                     
                    identical_seqs.gr<-GRanges(seqnames=identical_seqs.dt$first.seqnames, 
                    		   ranges=IRanges(start=identical_seqs.dt$first.start,
                            end=identical_seqs.dt$first.end),
                    		   strand="*")
                    
                    if(chrom != "all"){
                    	identical_seqs.gr<-subset(identical_seqs.gr, seqnames(identical_seqs.gr)==chrom)
                    }
                                
                    filtered_identical_seqs.gr <- identical_seqs.gr
                    
                    ####tile the genome
                    ##downloaded the .sizes file from https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/
                    
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
                    
                    tileHits<-findOverlaps(tiled_chr, filtered_identical_seqs.gr)
                    tileSplit<-split(subjectHits(tileHits), queryHits(tileHits))
                    
                    kurt_by_tile<-lapply(as.list(1:length(tiled_chr)) , function(x){
                    			     obj<-tileSplit[[as.character(x)]]
                    			     if(is.null(obj)){
                    				     return(0)
                    			     } else {
                    				     runs<-filtered_identical_seqs.gr[obj]
                    				     dist<-width(runs)
                    				     score<-(quantile(dist, 0.75) - (2 * quantile(dist, 0.50)) + quantile(dist, 0.25))/denominator
                    				     if(score=="NaN") score = 0
                    				     return(score)
                    			     }
                    			   })
                
                    tiled_chr$score<-unlist(kurt_by_tile)
                    
                    #drop windows that arent the correct width - the ends of chromosomes usually
                    
                    tiled_chr<-tiled_chr[width(tiled_chr)==window_size]
                    
                    export.wig(tiled_chr, paste0(outDir, species1, "_", names(species2)[i], "_", chrom, "_", window_size/1000, "kb_min_window.wig"))
                    
                    print("skewness calculated")
                }
          }
    }
}


#gc()