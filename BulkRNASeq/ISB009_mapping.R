library(data.table)
library(readxl)
library(stringr)

setwd('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB009')
files <- as.data.table(readxl::read_xlsx('/Volumes/ncatssctl/NGS_related/Sample Sheet/BulkRNA/ISB009_Library_prep_sample_sheet.xlsx', sheet=2))

files[,shortfile := tstrsplit(Sample, '_L')[1]]
files


dir <- '/data/NCATS_ifx/data/mRNASeq/ISB009'
samples <- unique(files$shortfile)
for (sample in samples){
  cat(paste0("java -jar $TRIMMOJAR PE -phred33 ", dir , '/', sample, '_L006_R1_001.fastq.gz ', dir, '/',
            sample, '_L006_R2_001.fastq.gz ', '-baseout ', sample, '.fastq.gz ',
             'ILLUMINACLIP:/usr/local/apps/trimmomatic/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:10 TRAILING:5 MAXINFO:50:0.97 MINLEN:36'), sep="\n")
  cat("\n")
}


## swarm -f trim_swarm_012919.sh -g 12 -t 8 --module trimmomatic #19325944
samples <- unique(files$SampleID)
for (sample in samples){
  cat(paste0('cd ', dir, ' && mkdir -p bam/' , 'Sample_' , sample , ' && STAR --runThreadN $SLURM_CPUS_PER_TASK --genomeDir /fdb/STAR_current/GENCODE/Gencode_human/release_27/genes-75 --sjdbOverhang 75 --outSAMunmapped Within --outFilterType BySJout --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --sjdbScore 1 --readFilesIn ', dir, '/' , sample , '_1P.fastq.gz ', dir, '/', sample , '_2P.fastq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix bam/Sample_' , sample ,'/' , sample, '_hg38'))
  cat("\n\n")
}
# swarm -f star_swarm_013119.sh -t 12 -g 40 --module STAR --time=8:00:00 




#

