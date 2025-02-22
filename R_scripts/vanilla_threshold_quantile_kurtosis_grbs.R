# Yout need to install R 3.4.2 in order to successfully compile this script because the updated version of GRanges objects don't support [[, as.list(), lapply(), or unlist() at the moment.
# source("https://bioconductor.org/biocLite.R")
# biocLite("CNEr")
# source("https://bioconductor.org/biocLite.R")
# biocLite("GenomicRanges")
library(GenomicRanges)
library(rtracklayer)
# install.packages('data.table')
library(data.table)
# install.packages('cpm')
library(cpm)
# install.packages('randtests')
library(randtests)
# install.packages('devtools')
# install.packages('remotes')
# library(remotes)
# devtools::install_github("ge11232002/CNEr")
# devtools::install_github('kloke/npsm')
# devtools::install_github("CshlSiepelLab/RPHAST")
library(npsm)
library(rphast)


##########Define variables

args<-commandArgs(TRUE)

species1<-'hg38'			## UCSC stle abbreviation (hg38, micMur3, etc)
species2<-list('calJac4'=30,'danRer10'=8,'equCab3'=12,'felCat9'=10,'galGal6'=8,'ornAna2'=7,'oryCun2'=10,'oviAri4'=11,'petMar3'=9,'susScr11'=10,'tarSyr2'=14,'thaSir1'=8,'triMan1'=10,'xenTro10'=10)	
window_list = list(as.numeric(15000), as.numeric(30000))    #list(as.numeric(15000), as.numeric(120000), as.numeric(30000), as.numeric(60000)) 	## in bp
ARL0=as.numeric(370)		## one of  370, 500, 600, 700, ..., 1000, 2000, 3000, ..., 10000, 20000, ..., 50000. see ?processStream
percentile=as.numeric(0.7)		## between 0-1, usually start with 0.7
#percent_list = list(as.numeric(100), as.numeric(98), as.numeric(96), as.numeric(96.6), as.numeric(90)) #it means the CNE identification threshold
percent_list = list(as.numeric(100))

#####################GRBS##########################

for(i in seq_along(species2))
{
    for(window_size in window_list)
    {
          for(percent in percent_list)
          {
                  # Import all the WIG files generated by the 2nd script named "binned_quantile_kurtosis_genome_ave.R"
                  #inDir<-paste0("/bucket/MillerU/Abrar/hg38-micMur3_outputs/", window_size/1000, "kb/binned_quantile_kurtosis_genome_ave", "_", percent, "pc_500col/")
                  inDir<-paste0("/bucket/MillerU/Abrar/hg38-", names(species2)[i], "_outputs/", window_size/1000, "kb/binned_quantile_kurtosis_genome_ave", "_", percent, "pc/")
                  
                  score_files<-dir(inDir, pattern=paste0(names(species2)[i], ".*", window_size/1000, "kb_min_window"), full.names=T)
                  
                  tiled_chr<-lapply(score_files, import.wig, genome=species1)
                  
                  tiled_chr_all<-Reduce(c, tiled_chr)
                  
                  edges<-lapply(tiled_chr, function(x){
                  		      fwrd.edges<-processStream(x$score, cpmType="Mann-Whitney", ARL0=ARL0)
                  		      rev.edges<-processStream(rev(x$score), cpmType="Mann-Whitney", ARL0=ARL0)
                  		      return(list("fwrd.edges"=fwrd.edges,
                  				  "rev.edges"=rev.edges))
                  })
                                 
                  boundaries<-mapply(x = tiled_chr,
                  		   y = edges, 
                  		   FUN = function(x,y){
                  			   x$fwrd.boundary=0
                  			   x$rev.boundary=0
                  			   x$fwrd.boundary[y$fwrd.edges$changePoints]=1
                  			   revx<-rev(x)
                  			   revx$rev.boundary[y$rev.edges$changePoints]=1
                  			   x$rev.boundary=rev(revx$rev.boundary)
                  			   x$score<-NULL
                  			   fwrd<-x
                  			   mcols(fwrd)<-NULL
                  			   fwrd$score<-x$fwrd.boundary
                  			   fwrd=subset(fwrd, score==1)
                  			   rev<-x
                  			   mcols(rev)<-NULL
                  			   rev$score<-x$rev.boundary
                  			   rev=subset(rev, score==1)
                  			   out<-list("fwrd"=fwrd, "rev"=rev)
                  			   return(out)
                  		   }, SIMPLIFY=FALSE)
                  
                  boundaries<-list("fwrd"=lapply(boundaries, function(x) x$fwrd),
                  		 "rev"=lapply(boundaries, function(x) x$rev))
                  
                  if(length(boundaries[[1]]) == 1){
                  	boundaries<-lapply(boundaries, function(x) split(x[[1]], seqnames(x[[1]])))
                  }
                  
                  ranges<-lapply(boundaries, function(x){
                  		       out<-Reduce(c, lapply(x, function(y){
                  						     if(length(y) < 2) return(GRanges())
                  						     ranges<-GRanges()
                  						     for(i in 1:length(y)){
                  							     if(i == length(y)) break
                  							     range<-GRanges(seqnames=seqnames(y)[i],
                  									    ranges=IRanges(start=start(y)[i],
                  											   end=end(y)[(i+1)]),
                  									    strand="*")
                  							     ranges<-c(ranges, range)
                  						     }
                  						     return(ranges)
                  				  }))
                  		       return(out)
                  		 })
                  
                  scores_by_ranges<-lapply(ranges, function(x){
                  				 out<-lapply(x, function(y){
                  						     scores<-subsetByOverlaps(tiled_chr_all, y)$score
                  						     return(scores)
                  })
                  				 return(out)
                  		 })
                  
                  mean_scores_by_ranges<-lapply(scores_by_ranges, function(x){
                  				      out<-do.call(c, lapply(x, mean))
                  		 })
                  
                  grbs_raw<-list()
                  
                  grbs_raw[["fwrd"]]<-reduce(ranges$fwrd[which(mean_scores_by_ranges$fwrd > quantile(mean_scores_by_ranges$fwrd, percentile))])
                  
                  grbs_raw[["rev"]]<-reduce(ranges$rev[which(mean_scores_by_ranges$rev > quantile(mean_scores_by_ranges$rev, percentile))])
                  
                  
                  combined_raw_grbs<-intersect(Reduce(c, grbs_raw$fwrd), Reduce(c, grbs_raw$rev))
    
                  
                  #outDir<-paste0("/bucket/MillerU/Abrar/hg38-micMur3_outputs/", window_size/1000, "kb/binned_quantile_kurtosis_genome_ave_grbs", "_", percent, "pc_500col/")
                  outDir<-paste0("/bucket/MillerU/Abrar/hg38-", names(species2)[i], "_outputs/", window_size/1000, "kb/binned_quantile_kurtosis_genome_ave_grbs", "_", percent, "pc/")
                  
                  dir.create(outDir, showWarnings=F)
                  
                  export.bed(combined_raw_grbs, paste0(outDir, species1, "_", names(species2)[i], "_binned_kurtosis_arl_", ARL0, "_", percentile, "_pc_", window_size/1000, "_kb_windows_grbs.bed"))
                  
                  message("Script Complete")
                  #gc()
    
          }
    }
}

