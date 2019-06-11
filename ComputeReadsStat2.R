#Eric Fournier 2019-03-27
#Script appelé par MiSeqStat5.py pour calculer les metrics d'une run MiSeq

library("ShortRead")
library("tidyverse")


#pour ligne de commande
args <- commandArgs(TRUE)

#le repertoire contenant les fastq
fastqpath <- args[1]

#repertoire de calcul temporaire
tempwd <- args[2]

fastqpath <- file.path(fastqpath)

setwd(tempwd)

outfile <- file.path(getwd(),"fastqStat.txt")
#print(outfile)

dummy_fastq_file <- system.file(package = "ShortRead","extdata", "E-MTAB-1147","ERR127302_1_subset.fastq.gz")
dummy_fastq_object <- readFastq(dummy_fastq_file)
#print(dummy_fastq_object)

qual_index <- encoding(quality(dummy_fastq_object))
#print(qual_index)

fastq_dir <- dir(fastqpath,"*fastq.gz", full = TRUE)

qaSummary <- qa(fastq_dir, type="fastq")
#print(qaSummary)

readCount_df <- qaSummary[["readCounts"]]
#print(readCount_df)
readCount_df <- cbind(readCount_df, row.names(readCount_df))
colnames(readCount_df)[4] <- "readname"
#print(readCount_df)

baseQual_df <- qaSummary[["baseQuality"]]
#print(baseQual_df)
baseQual_df <- mutate(baseQual_df, IntQual=qual_index[as.character(score)])
#print(baseQual_df)

fastq_name_list <- distinct(baseQual_df, lane)
fastq_name_list <- fastq_name_list[,1]
#print(fastq_name_list)

#somme des nucleotides a travers tous les fastq
sum_nt_all_q_allspec <- 0
#somme des nucleotides avec minimum q30 a travers tous les fastq
sum_nt_min_q30_allspec <- 0
#somme des reads a travers tous les fastq
sum_reads_allspec <- 0

for(fastqname in fastq_name_list){
  
  nb_reads <- filter(readCount_df, readname == fastqname)
  nb_reads <- nb_reads[,1]
  sum_reads_allspec <- sum_reads_allspec + nb_reads
  
  #toutes les lignes qualite pour ce fastq
  wholeQ_df <- filter(baseQual_df, lane == fastqname)
  #vector avec les comptes de nucleotides pour toutes les valeurs de qualité
  wholeQ_nt_count <- wholeQ_df[,2]
  #somme des nucleotides
  wholeQ_nt_count_sum <- sum(wholeQ_nt_count)
  #on incremente la somme des nucleotides totales pour cette run
  sum_nt_all_q_allspec <- sum_nt_all_q_allspec + wholeQ_nt_count_sum
  
  #idem a ci-dessus mais pour les nucleotides ayant une valeur Q minimum de 30
  minQ30_df <- filter(baseQual_df, lane == fastqname & IntQual > 29)
  minQ30_nt_count <- minQ30_df[,2]
  minQ30_nt_count_sum <- sum(minQ30_nt_count)
  sum_nt_min_q30_allspec <- sum_nt_min_q30_allspec + minQ30_nt_count_sum
  
  short_fastqname <- str_replace(fastqname,"(\\S+)_\\S+_\\S+_(\\S+)_\\S+.fastq.gz","\\1_\\2")
  
  #proportion de nucleotide ayant une valeur Q de minimum 30 pour ce fastq
  perc_nt_min_q30 <- round((minQ30_nt_count_sum/wholeQ_nt_count_sum)*100,0)
  #print(perc_nt_min_q30)
  
  #On enregistre les metrics pour ce specimen
  write(paste0(short_fastqname,"\t",perc_nt_min_q30,"\t",nb_reads,"\t",wholeQ_nt_count_sum),file=outfile,append = TRUE)
  
}

#Metrics de la run
perc_nt_min_q30_allspec <- round((sum_nt_min_q30_allspec/sum_nt_all_q_allspec)*100,0)
write(paste0("RUN","\t",perc_nt_min_q30_allspec,"\t",sum_reads_allspec,"\t",sum_nt_min_q30_allspec),file=outfile,append = TRUE)







