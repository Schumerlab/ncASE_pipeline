#usage: replicate_1_par_1 replicate_1_par_2 ... filter read_length output
#Run through ncASE_pipeline, launched from ncASE_parallel.pl:
#Rscript ncASE_cmd_I-2020.R $ase_string $allow_zero_counts $bias_threshold $read_length $summary_file
#
#Assumptions: 
#  - input files sorted by transcript and position
#  - vcf info from two alignments has been reorganized s.t.
#    the reciprocal ref and alt alleles are made the same 
#vi-2020: compatible with MS ncASE_pipeline

library(tidyverse)

args<-commandArgs(TRUE)

#We assume the input will be in sorted order (of both transcript and position)
snp_weight <- function(x, read_length) {
   dist_mat <- 1 - abs(outer(x$position, x$position, "-"))/read_length
   dist_mat[dist_mat < 0] <- 0
   weights <- 1/colSums(dist_mat)
   return(data.frame(transcript=x$transcript, position=x$position, weight=weights, stringsAsFactors=FALSE))
}

if(length(args) < 6)
   stop("Need at least one pair of replicate files and output filename")
if((length(args)-1) %% 2 != 1)
   stop("Need two files per replicate")

allow_zero_counts <- as.numeric(args[length(args)-3])

read_length <- as.numeric(args[length(args)-1]) 

max_replicate <- (length(args)-4)/2

transcript_count_list <- vector("list", max_replicate)

for(replicate in 1:max_replicate) {
   print(replicate)

   cat("Processing replicate", replicate, "\n")

   par1_cols <- c("transcript", "position", "par1ref", "par1alt", "par1ref_count", "par1alt_count")
   par2_cols <- c("transcript", "position", "par2ref", "par2alt", "par2ref_count", "par2alt_count")
   par_col_classes <- c("character", "integer", "character", "character", "integer", "integer")
   
   par1 <- read.table(args[2*replicate-1], quote="\"", col.names=par1_cols, colClasses=par_col_classes)
   cat(paste(nrow(par1), "par1 SNPs from input file\n"))
   par2 <- read.table(args[2*replicate], quote="\"", col.names=par2_cols, colClasses=par_col_classes)
   cat(paste(nrow(par2), "par2 SNPs from input file\n"))

   # If allow_zero_counts, filter for par[12]ref_count > 0 
   # && par[12]alt_count > 0, and set up a column as a join key, 
   # eliminating the redundant transcript and position columns:
   if(!allow_zero_counts) {
     par1 <- par1 %>% filter(par1ref_count * par1alt_count > 0) %>%
        mutate(variant=paste(transcript, position, sep="_")) %>%
        select(-c(transcript, position))
     cat(paste(nrow(par1), "par1 SNPs after non-zero ref*alt count filter\n"))
     par2 <- par2 %>% filter(par2ref_count * par2alt_count > 0) %>%
        mutate(variant=paste(transcript, position, sep="_")) %>%
        select(-c(transcript, position))
     cat(paste(nrow(par2), "par2 SNPs after non-zero ref*alt count filter\n"))
   } else {
     par1 <- par1 %>%
             mutate(variant=paste(transcript, position, sep="_")) %>%
             select(-c(transcript, position))
     par2 <- par2 %>% 
             mutate(variant=paste(transcript, position, sep="_")) %>%
             select(-c(transcript, position))
     cat("Note: Zero counts retained.\n")
   }

   #Perform an inner join on the sets of SNPs identified by mapping to each
   # parent:
   snps <- inner_join(par1, par2, by="variant")
   cat(paste(nrow(snps), "SNPs after inner join of the two parent mappings\n"))

   #Filter for SNPs with matching ref and alt on the two parental
   # transcriptomes (e.g. par1ref == par2ref && par1alt == par2alt:
   #Note: pipeline already does this, serves as a safety check 
   snps <- snps %>% filter(par1ref == par2ref & par1alt == par2alt) %>%
      separate(variant, into=c("transcript", "position"), sep="_", remove=FALSE) %>%
      mutate(position=as.integer(position)) %>%
      arrange(transcript, position)
   cat(paste(nrow(snps), "SNPs after reciprocal parental allele filter\n"))

   #Set thresholds for proximity of reciprocal mapping counts:
   filter_threshold <- as.numeric(args[length(args)-2])
   #The old threshold+filter combination is equivalent to the following:
   #low_threshold <- (1 - filter_threshold)/(1 + filter_threshold)
   #high_threshold <- (1 + filter_threshold)/(1 - filter_threshold)
   #New thresholds, representing % absolute deviance in either direction:
   low_threshold <- 1 - filter_threshold
   high_threshold <- 1 + filter_threshold
   #Old thresholds for old filters:
   #low_threshold <- 0.5 - filter_threshold/2
   #high_threshold <- 0.5 + filter_threshold/2

   # ncASE_pipeline _ASE_count input has the same allele in the ref
   # position for both par1 and par2, regardless of alignment genome
   # so, par1ref and par2ref should be compared, and par1alt and par2alt
   # should be compared
   snps <- snps %>% filter(par1ref_count > low_threshold * par2ref_count)
   cat(paste(nrow(snps), "SNPs after high-pass par1ref filter\n"))
   snps <- snps %>% filter(par1ref_count < high_threshold * par2ref_count)
   cat(paste(nrow(snps), "SNPs after adding low-pass par1ref filter\n"))
   snps <- snps %>% filter(par2alt_count > low_threshold * par1alt_count)
   cat(paste(nrow(snps), "SNPs after adding high-pass par1alt filter\n"))
   snps <- snps %>% filter(par2alt_count < high_threshold * par1alt_count)
   cat(paste(nrow(snps), "SNPs after all reciprocal mapping counts filters\n"))
   print(snps)

   # Calculate a per-SNP weight based on distance in the transcript:
   transcripts <- unique(snps$transcript)
   weight_df <- data.frame(transcript=c(), position=c(), weight=c())
   for (tx in transcripts) {
      weight_df <- rbind(weight_df, snp_weight(snps[snps$transcript == tx,], read_length))
   }
   cat(paste("Calculated SNP weights for", length(transcripts), "transcripts\n"))
   #cat(paste("snps$transcript class", class(snps$transcript), "snps$position class", class(snps$position), "weight_df$transcript class", class(weight_df$transcript), "weight_df$position class", class(weight_df$position), "\n"))
   snps <- full_join(snps, weight_df, by=c("transcript", "position"))
   cat(paste("Joined SNP weights into main data.frame for replicate", replicate, "\n"))

   # Calculate weighted counts for each transcript, get 
   # average of ref and alt counts from each par[12]_ASE_counts
   transcript_counts <- snps %>% group_by(transcript) %>%
      summarize(par1ref_counts=sum(par1ref_count*weight),
         par1alt_counts=sum(par1alt_count*weight),
         par2ref_counts=sum(par2ref_count*weight),
         par2alt_counts=sum(par2alt_count*weight)) %>%
      transmute(transcript=transcript,
         par1_counts=(par1ref_counts+par2ref_counts)/2,
         par2_counts=(par1alt_counts+par2alt_counts)/2)
   cat(paste("Summarized per-parent transcript counts for", nrow(transcript_counts), "transcripts for replicate", replicate, "\n"))

   print(transcript_counts)
   colnames(transcript_counts)[2:3] <- c(paste("rep", replicate, "par1", "counts", sep="_"), paste("rep", replicate, "par2", "counts", sep="_"))
   cat("Renamed count columns to include replicate\n")
   transcript_count_list[[replicate]] <- transcript_counts
}

# Join all the transcript count tables together and 
# output in summary file
i <- 1
final_transcript_counts <- transcript_count_list[[i]]
cat(paste(nrow(final_transcript_counts), "transcripts from replicate", i, "\n"))
i <- i + 1
while (i <= max_replicate) {
   final_transcript_counts <- final_transcript_counts %>% full_join(transcript_count_list[[i]], by="transcript")
   cat(paste(nrow(final_transcript_counts), "transcripts after joining with replicate", i, "\n"))
   i <- i + 1
}

write.table(final_transcript_counts, file=args[length(args)], quote=FALSE, sep="\t", row.names=FALSE)
