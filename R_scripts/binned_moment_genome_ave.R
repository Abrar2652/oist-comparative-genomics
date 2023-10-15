library(GenomicRanges)
library(rtracklayer)
library(data.table)
library(npsm)
library(rphast)
#install.packages('e1071')
library(e1071)   


##################################################################################################

args<-commandArgs(TRUE)

window_list = list(as.numeric(30000), as.numeric(15000), as.numeric(120000), as.numeric(60000))  	##in bp
species1<-'hg19'			## UCSC stle abbreviation (hg19, canFam3, etc)	
species2<-'mm10'		
## either "all" or UCSC style notation (chr1, chr2, chrX, etc)
chrom_list=list('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY')
#chrom_list='all'
percent_list = list(as.numeric(100)) #it means the CNE identification threshold

for(window_size in window_list)
{        
    for(percent in percent_list)
    {
        # Only for Gorilla
        #outDir1<-paste0("/bucket/MillerU/Abrar/hg19-mm10_outputs/", window_size/1000, "kb/binned_mean_genome_ave","_",percent,"pc_500col/")
        #outDir2<-paste0("/bucket/MillerU/Abrar/hg19-mm10_outputs/", window_size/1000, "kb/binned_variance_genome_ave","_",percent,"pc_500col/")
        #outDir3<-paste0("/bucket/MillerU/Abrar/hg19-mm10_outputs/", window_size/1000, "kb/binned_skewness_genome_ave","_",percent,"pc_500col/")
        #outDir4<-paste0("/bucket/MillerU/Abrar/hg19-mm10_outputs/", window_size/1000, "kb/binned_excess_kurtosis_genome_ave","_",percent,"pc_500col/")
        #outDir5<-paste0("/bucket/MillerU/Abrar/hg19-mm10_outputs/", window_size/1000, "kb/binned_sixth_moment_genome_ave","_",percent,"pc_500col/")
        # Others
        outDir1<-paste0("/bucket/MillerU/Abrar/hg19-mm10_outputs/", window_size/1000, "kb/binned_mean_genome_ave","_",percent,"pc/")
        outDir2<-paste0("/bucket/MillerU/Abrar/hg19-mm10_outputs/", window_size/1000, "kb/binned_variance_genome_ave","_",percent,"pc/")
        outDir3<-paste0("/bucket/MillerU/Abrar/hg19-mm10_outputs/", window_size/1000, "kb/binned_skewness_genome_ave","_",percent,"pc/")
        outDir4<-paste0("/bucket/MillerU/Abrar/hg19-mm10_outputs/", window_size/1000, "kb/binned_excess_kurtosis_genome_ave","_",percent,"pc/")
        outDir5<-paste0("/bucket/MillerU/Abrar/hg19-mm10_outputs/", window_size/1000, "kb/binned_sixth_moment_genome_ave","_",percent,"pc/")
        
        if(!dir.exists(file.path(outDir1))){
          dir.create(file.path(outDir1))
        }
        if(!dir.exists(file.path(outDir2))){
          dir.create(file.path(outDir2))
        }
        if(!dir.exists(file.path(outDir3))){
          dir.create(file.path(outDir3))
        }
        if(!dir.exists(file.path(outDir4))){
          dir.create(file.path(outDir4))
        }
        if(!dir.exists(file.path(outDir5))){
          dir.create(file.path(outDir5))
        }  
          
        for (chrom in chrom_list)
        {    
            ## read in the identical sequence locations (tsv files) - output of filter_exons_repeats.py
            ## get_identical_seq_loc_hg19-mm10     
            identical_seqs_fn<-paste0("/bucket/MillerU/Abrar/hg19-mm10_outputs/filtered_pcs_hg19-mm10/", chrom)
            
            if(chrom == "all"){
              identical_seqs_files<-dir(identical_seqs_fn, pattern=paste0(species1, "_", species2, ".*all*"), full.names=T)
            } else {
              identical_seqs_files<-dir(identical_seqs_fn, pattern=paste0(species1, "_", species2, ".*chr*"), full.names=T)
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
            ##downloaded the .sizes file from https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/
            
            genome_info<-read.table(paste0("/bucket/MillerU/Abrar/bigZips/", species1, ".chrom.sizes"))
            
            genome_seqlengths<-genome_info$V2
            names(genome_seqlengths)<-genome_info$V1
            
            tiled_genome<-tileGenome(genome_seqlengths, tilewidth=window_size, cut.last.tile.in.chrom=T)
            
            if(chrom=="all"){
               tiled_chr<-tiled_genome
             	 tiled_chr1<-tiled_genome
               tiled_chr2<-tiled_genome
               tiled_chr3<-tiled_genome
               tiled_chr4<-tiled_genome
               tiled_chr5<-tiled_genome
            } else {
               tiled_chr<-tiled_genome[seqnames(tiled_genome)==chrom]
           	   tiled_chr1<-tiled_genome[seqnames(tiled_genome)==chrom]
               tiled_chr2<-tiled_genome[seqnames(tiled_genome)==chrom]
               tiled_chr3<-tiled_genome[seqnames(tiled_genome)==chrom]
               tiled_chr4<-tiled_genome[seqnames(tiled_genome)==chrom]
               tiled_chr5<-tiled_genome[seqnames(tiled_genome)==chrom]
            }
            
            ####get distribution in each window and calculate kurtosis
            
            tileHits<-findOverlaps(tiled_chr, filtered_identical_seqs.gr)
            tileSplit<-split(subjectHits(tileHits), queryHits(tileHits))
            
            # mean
            kurt_by_tile1<-lapply(as.list(1:length(tiled_chr1)) , function(x){
              obj<-tileSplit[[as.character(x)]]
              if(is.null(obj)){
                return(0)
              } else {
                runs<-filtered_identical_seqs.gr[obj]
                dist<-width(runs)
                score<-mean(dist)
                if(score=="NaN") score = 0
                return(score)
              }
            })
            
            # variance
            kurt_by_tile2<-lapply(as.list(1:length(tiled_chr2)) , function(x){
              obj<-tileSplit[[as.character(x)]]
              if(is.null(obj)){
                return(0)
              } else {
                runs<-filtered_identical_seqs.gr[obj]
                dist<-width(runs)
                score<-var(dist) 
                if(is.na(score)) score = 0
                return(score)
              }
            })
            
            # skewness
            kurt_by_tile3<-lapply(as.list(1:length(tiled_chr3)) , function(x){
              obj<-tileSplit[[as.character(x)]]
              if(is.null(obj)){
                return(0)
              } else {
                runs<-filtered_identical_seqs.gr[obj]
                dist<-width(runs)
                score<-skewness(dist)
                if(score=="NaN") score = 0
                return(score)
              }
            })
            
            # excess kurtosis or Fisher's Kurtosis
            kurt_by_tile4<-lapply(as.list(1:length(tiled_chr4)) , function(x){
              obj<-tileSplit[[as.character(x)]]
              if(is.null(obj)){
                return(0)
              } else {
                runs<-filtered_identical_seqs.gr[obj]
                dist<-width(runs)
                score<-kurtosis(dist)
                if(score=="NaN") score = 0
                return(score)
              }
            })
            
            # 6th central moment
            kurt_by_tile5<-lapply(as.list(1:length(tiled_chr5)) , function(x){
              obj<-tileSplit[[as.character(x)]]
              if(is.null(obj)){
                return(0)
              } else {
                runs<-filtered_identical_seqs.gr[obj]
                dist<-width(runs)
                score<-moment(dist, order=6, center=TRUE) 
                if(score=="NaN") score = 0
                return(score)
              }
            })
            				
            
            tiled_chr1$score<-unlist(kurt_by_tile1)
            tiled_chr2$score<-unlist(kurt_by_tile2)
            tiled_chr3$score<-unlist(kurt_by_tile3)
            tiled_chr4$score<-unlist(kurt_by_tile4)
            tiled_chr5$score<-unlist(kurt_by_tile5)
            
            #drop windows that arent the correct width - the ends of chromosomes usually
            
            tiled_chr1<-tiled_chr1[width(tiled_chr1)==window_size]
            tiled_chr2<-tiled_chr2[width(tiled_chr2)==window_size]
            tiled_chr3<-tiled_chr3[width(tiled_chr3)==window_size]
            tiled_chr4<-tiled_chr4[width(tiled_chr4)==window_size]
            tiled_chr5<-tiled_chr5[width(tiled_chr5)==window_size]
            
            export.wig(tiled_chr1, paste0(outDir1, species1, "_", species2, "_", chrom, "_", window_size/1000, "kb_min_window.wig"))
            export.wig(tiled_chr2, paste0(outDir2, species1, "_", species2, "_", chrom, "_", window_size/1000, "kb_min_window.wig"))
            export.wig(tiled_chr3, paste0(outDir3, species1, "_", species2, "_", chrom, "_", window_size/1000, "kb_min_window.wig"))
            export.wig(tiled_chr4, paste0(outDir4, species1, "_", species2, "_", chrom, "_", window_size/1000, "kb_min_window.wig"))
            export.wig(tiled_chr5, paste0(outDir5, species1, "_", species2, "_", chrom, "_", window_size/1000, "kb_min_window.wig"))
        
            print("moments calculated")
        }
    }  
}










