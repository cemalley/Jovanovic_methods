library(data.table)
library(readxl)
library(stringr)

setwd('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB025/')
files <- as.data.table(readxl::read_xlsx('ISB025_lane_concatenating.xlsx', sheet=1))

samples <- unique(files$Sample)

# ~~~~~~~PC_siRNA~~~~~~~~-------------------
setwd('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB025/')
files <- as.data.table(readxl::read_xlsx('ISB025_lane_concatenating.xlsx'))

samples <- unique(files$Sample)

# lane concatenation-----
for (sample in samples){
  cat(paste0('mkdir ', sample), sep='\n')
}


for (sample in samples){
  rows <- grep(sample, files$Sample)
  cat(sample)
  cat('\t')
  cat(rows)
  cat('\n')
}

rows <- readxl::read_xlsx('ISB025_lane_concatenating.xlsx', sheet=2)
rows <- as.data.table(rows)
rows

samples <- unique(rows$Sample)

#R1------
for (sample in samples){
  files.row <- grep(sample, rows$Sample)
  files.temp <- unlist(as.numeric(str_split_fixed(rows$Rows[files.row], " ", Inf)))
  
  files.temp.out <- c()
  cat('cat ', sep='') 
  
  for (row in files.temp){
    files.temp[row] <- files$File[row]
    dir.row <- grep(files.temp[row],files$File)
    files.temp.out[row] <- paste0('/data/NCATS_ifx/data/mRNASeq/ISB025/ISB025_PC_siRNA/Raw_data/FASTQ/' ,files$Folder[dir.row], '/', files$File[dir.row])
    cat(paste0(files.temp.out[row], '_R1_001.fastq.gz '))
  }
  
  out.sample <- files$Sample[dir.row]
  
  cat(paste0('> ', '/data/NCATS_ifx/data/mRNASeq/ISB025/ISB025_PC_siRNA/Raw_data/concatenated_lanes/', out.sample, '/', out.sample, '_R1_001.fastq.gz'))
  cat("\n\n", sep="")
  
  rm(files.temp.out, files.row, files.temp, dir.row)
}
# R2----
for (sample in samples){
  files.row <- grep(sample, rows$Sample)
  files.temp <- unlist(as.numeric(str_split_fixed(rows$Rows[files.row], " ", Inf)))
  
  files.temp.out <- c()
  cat('cat ', sep='') 
  
  for (row in files.temp){
    files.temp[row] <- files$File[row]
    dir.row <- grep(files.temp[row],files$File)
    files.temp.out[row] <- paste0('/data/NCATS_ifx/data/mRNASeq/ISB025/ISB025_PC_siRNA/Raw_data/FASTQ/' ,files$Folder[dir.row], '/', files$File[dir.row])
    cat(paste0(files.temp.out[row], '_R2_001.fastq.gz '))
  }
  
  out.sample <- files$Sample[dir.row]
  
  cat(paste0('> ', '/data/NCATS_ifx/data/mRNASeq/ISB025/ISB025_PC_siRNA/Raw_data/concatenated_lanes/', out.sample, '/', out.sample, '_R2_001.fastq.gz'))
  cat("\n\n", sep="")
  
  rm(files.temp.out, files.row, files.temp, dir.row)
}
# swarm -f fastq_cat_081820.sh -t 1 -g 1 --time=8:00:00 #

# mapping-------
samples <- unique(files$Sample)

dir <- '/data/NCATS_ifx/data/mRNASeq/ISB025/ISB025_PC_siRNA'
for (sample in samples){
  cat(paste0('cd ',dir,'/trimmomatic && java -jar $TRIMMOJAR PE -phred33 ', dir ,'/Raw_data/concatenated_lanes/',sample,'/', sample,'_R1_001.fastq.gz ', dir, '/Raw_data/concatenated_lanes/',sample, '/',
             sample, '_R2_001.fastq.gz ', '-baseout ', sample, '.fastq.gz ',
             'ILLUMINACLIP:/usr/local/apps/trimmomatic/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:10 TRAILING:5 MAXINFO:50:0.97 MINLEN:36 && cd ',dir,'/BAM && STAR --runThreadN $SLURM_CPUS_PER_TASK --genomeDir /fdb/STAR_current/GENCODE/Gencode_human/release_27/genes-150 --sjdbOverhang 150 --outSAMunmapped Within --outFilterType BySJout --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --sjdbScore 1 --readFilesIn ',dir,'/trimmomatic/',sample,'_1P.fastq.gz ',dir,'/trimmomatic/',sample,'_2P.fastq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ', sample,'_hg38 && htseq-count -f bam -r pos -s no -t exon -m union ',dir,'/BAM/',sample,'_hg38Aligned.sortedByCoord.out.bam /fdb/GENCODE/Gencode_human/release_27/gencode.v27.annotation.gtf --idattr gene_name > ',dir,'/htseq/',sample,'_htseq_counts.txt'), sep="\n")
  cat("\n")
}

# swarm -g 40 -t 12 --time=48:00:00 --module=trimmomatic,STAR,htseq -f ISB025_PC_siRNA_pipeline.sh #



# ~~~~~~~PC_LA_WNT~~~~~~~~-------------------
setwd('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB025/')
files <- as.data.table(readxl::read_xlsx('ISB025_lane_concatenating.xlsx', sheet=3))

samples <- unique(files$Sample)

# lane concatenation-----
for (sample in samples){
  cat(paste0('mkdir ', sample), sep='\n')
}


for (sample in samples){
  rows <- grep(sample, files$Sample)
  cat(sample)
  cat('\t')
  cat(rows)
  cat('\n')
}

rows <- readxl::read_xlsx('ISB025_lane_concatenating.xlsx', sheet=4)
rows <- as.data.table(rows)
rows

samples <- unique(rows$Sample)

#R1------
for (sample in samples){
  files.row <- grep(sample, rows$Sample)
  files.temp <- unlist(as.numeric(str_split_fixed(rows$Rows[files.row], " ", Inf)))
  
  files.temp.out <- c()
  cat('cat ', sep='') 
  
  for (row in files.temp){
    files.temp[row] <- files$File[row]
    dir.row <- grep(files.temp[row],files$File)
    files.temp.out[row] <- paste0('/data/NCATS_ifx/data/mRNASeq/ISB025/ISB025_PC_LA_WNT/Raw_data/FASTQ/' ,files$Folder[dir.row], '/', files$File[dir.row])
    cat(paste0(files.temp.out[row], '_R1_001.fastq.gz '))
  }
  
  out.sample <- files$Sample[dir.row]
  
  cat(paste0('> ', '/data/NCATS_ifx/data/mRNASeq/ISB025/ISB025_PC_LA_WNT/Raw_data/concatenated_lanes/', out.sample, '/', out.sample, '_R1_001.fastq.gz'))
  cat("\n\n", sep="")
  
  rm(files.temp.out, files.row, files.temp, dir.row)
}
# R2----
for (sample in samples){
  files.row <- grep(sample, rows$Sample)
  files.temp <- unlist(as.numeric(str_split_fixed(rows$Rows[files.row], " ", Inf)))
  
  files.temp.out <- c()
  cat('cat ', sep='') 
  
  for (row in files.temp){
    files.temp[row] <- files$File[row]
    dir.row <- grep(files.temp[row],files$File)
    files.temp.out[row] <- paste0('/data/NCATS_ifx/data/mRNASeq/ISB025/ISB025_PC_LA_WNT/Raw_data/FASTQ/' ,files$Folder[dir.row], '/', files$File[dir.row])
    cat(paste0(files.temp.out[row], '_R2_001.fastq.gz '))
  }
  
  out.sample <- files$Sample[dir.row]
  
  cat(paste0('> ', '/data/NCATS_ifx/data/mRNASeq/ISB025/ISB025_PC_LA_WNT/Raw_data/concatenated_lanes/', out.sample, '/', out.sample, '_R2_001.fastq.gz'))
  cat("\n\n", sep="")
  
  rm(files.temp.out, files.row, files.temp, dir.row)
}
# swarm -f fastq_cat_081820.sh -t 1 -g 1 --time=8:00:00 #

# mapping-------
samples <- unique(files$Sample)

dir <- '/data/NCATS_ifx/data/mRNASeq/ISB025/ISB025_PC_LA_WNT'
for (sample in samples){
  cat(paste0('cd ',dir,'/trimmomatic && java -jar $TRIMMOJAR PE -phred33 ', dir ,'/Raw_data/concatenated_lanes/',sample,'/', sample,'_R1_001.fastq.gz ', dir, '/Raw_data/concatenated_lanes/',sample, '/',
             sample, '_R2_001.fastq.gz ', '-baseout ', sample, '.fastq.gz ',
             'ILLUMINACLIP:/usr/local/apps/trimmomatic/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:10 TRAILING:5 MAXINFO:50:0.97 MINLEN:36 && cd ',dir,'/BAM && STAR --runThreadN $SLURM_CPUS_PER_TASK --genomeDir /fdb/STAR_current/GENCODE/Gencode_human/release_27/genes-150 --sjdbOverhang 150 --outSAMunmapped Within --outFilterType BySJout --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --sjdbScore 1 --readFilesIn ',dir,'/trimmomatic/',sample,'_1P.fastq.gz ',dir,'/trimmomatic/',sample,'_2P.fastq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ', sample,'_hg38 && htseq-count -f bam -r pos -s no -t exon -m union ',dir,'/BAM/',sample,'_hg38Aligned.sortedByCoord.out.bam /fdb/GENCODE/Gencode_human/release_27/gencode.v27.annotation.gtf --idattr gene_name > ',dir,'/htseq/',sample,'_htseq_counts.txt'), sep="\n")
  cat("\n")
}

# swarm -g 40 -t 12 --time=48:00:00 --module=trimmomatic,STAR,htseq -f ISB025_PC_LA_WNT_pipeline.sh #

# ~~~~~~~Vukasin~~~~~~~~-------------------
setwd('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB025/')
files <- as.data.table(readxl::read_xlsx('ISB025_lane_concatenating.xlsx', sheet=5))

samples <- unique(files$Sample)

# lane concatenation-----
for (sample in samples){
  cat(paste0('mkdir ', sample), sep='\n')
}


for (sample in samples){
  rows <- grep(sample, files$Sample)
  cat(sample)
  cat('\t')
  cat(rows)
  cat('\n')
}

rows <- readxl::read_xlsx('ISB025_lane_concatenating.xlsx', sheet=6)
rows <- as.data.table(rows)
rows

samples <- unique(rows$Sample)

#R1------
for (sample in samples){
  files.row <- grep(sample, rows$Sample)
  files.temp <- unlist(as.numeric(str_split_fixed(rows$Rows[files.row], " ", Inf)))
  
  files.temp.out <- c()
  cat('cat ', sep='') 
  
  for (row in files.temp){
    files.temp[row] <- files$File[row]
    dir.row <- grep(files.temp[row],files$File)
    files.temp.out[row] <- paste0('/data/NCATS_ifx/data/mRNASeq/ISB025/ISB025_Vukasin/Raw_data/FASTQ/' ,files$Folder[dir.row], '/', files$File[dir.row])
    cat(paste0(files.temp.out[row], '_R1_001.fastq.gz '))
  }
  
  out.sample <- files$Sample[dir.row]
  
  cat(paste0('> ', '/data/NCATS_ifx/data/mRNASeq/ISB025/ISB025_Vukasin/Raw_data/concatenated_lanes/', out.sample, '/', out.sample, '_R1_001.fastq.gz'))
  cat("\n\n", sep="")
  
  rm(files.temp.out, files.row, files.temp, dir.row)
}
# R2----
for (sample in samples){
  files.row <- grep(sample, rows$Sample)
  files.temp <- unlist(as.numeric(str_split_fixed(rows$Rows[files.row], " ", Inf)))
  
  files.temp.out <- c()
  cat('cat ', sep='') 
  
  for (row in files.temp){
    files.temp[row] <- files$File[row]
    dir.row <- grep(files.temp[row],files$File)
    files.temp.out[row] <- paste0('/data/NCATS_ifx/data/mRNASeq/ISB025/ISB025_Vukasin/Raw_data/FASTQ/' ,files$Folder[dir.row], '/', files$File[dir.row])
    cat(paste0(files.temp.out[row], '_R2_001.fastq.gz '))
  }
  
  out.sample <- files$Sample[dir.row]
  
  cat(paste0('> ', '/data/NCATS_ifx/data/mRNASeq/ISB025/ISB025_Vukasin/Raw_data/concatenated_lanes/', out.sample, '/', out.sample, '_R2_001.fastq.gz'))
  cat("\n\n", sep="")
  
  rm(files.temp.out, files.row, files.temp, dir.row)
}
# swarm -f fastq_cat_081820.sh -t 1 -g 1 --time=8:00:00 #

# mapping-------
samples <- unique(files$Sample)

dir <- '/data/NCATS_ifx/data/mRNASeq/ISB025/ISB025_Vukasin'
for (sample in samples){
  cat(paste0('cd ',dir,'/trimmomatic && java -jar $TRIMMOJAR PE -phred33 ', dir ,'/Raw_data/concatenated_lanes/',sample,'/', sample,'_R1_001.fastq.gz ', dir, '/Raw_data/concatenated_lanes/',sample, '/',
             sample, '_R2_001.fastq.gz ', '-baseout ', sample, '.fastq.gz ',
             'ILLUMINACLIP:/usr/local/apps/trimmomatic/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:10 TRAILING:5 MAXINFO:50:0.97 MINLEN:36 && cd ',dir,'/BAM && STAR --runThreadN $SLURM_CPUS_PER_TASK --genomeDir /fdb/STAR_current/GENCODE/Gencode_human/release_27/genes-150 --sjdbOverhang 150 --outSAMunmapped Within --outFilterType BySJout --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --sjdbScore 1 --readFilesIn ',dir,'/trimmomatic/',sample,'_1P.fastq.gz ',dir,'/trimmomatic/',sample,'_2P.fastq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ', sample,'_hg38 && htseq-count -f bam -r pos -s no -t exon -m union ',dir,'/BAM/',sample,'_hg38Aligned.sortedByCoord.out.bam /fdb/GENCODE/Gencode_human/release_27/gencode.v27.annotation.gtf --idattr gene_name > ',dir,'/htseq/',sample,'_htseq_counts.txt'), sep="\n")
  cat("\n")
}

# swarm -g 40 -t 12 --time=48:00:00 --module=trimmomatic,STAR,htseq -f ISB025_Vukasin_pipeline.sh #


# ~~~~~~~CT~~~~~~~~-------------------
setwd('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB025/')
files <- as.data.table(readxl::read_xlsx('ISB025_lane_concatenating.xlsx', sheet=7))

samples <- unique(files$Sample)

# lane concatenation-----
for (sample in samples){
  cat(paste0('mkdir ', sample), sep='\n')
}


for (sample in samples){
  rows <- grep(sample, files$Sample)
  cat(sample)
  cat('\t')
  cat(rows)
  cat('\n')
}

rows <- readxl::read_xlsx('ISB025_lane_concatenating.xlsx', sheet=8)
rows <- as.data.table(rows)
rows

samples <- unique(rows$Sample)

#R1------
for (sample in samples){
  files.row <- grep(sample, rows$Sample)
  files.temp <- unlist(as.numeric(str_split_fixed(rows$Rows[files.row], " ", Inf)))
  
  files.temp.out <- c()
  cat('cat ', sep='') 
  
  for (row in files.temp){
    files.temp[row] <- files$File[row]
    dir.row <- grep(files.temp[row],files$File)
    files.temp.out[row] <- paste0('/data/NCATS_ifx/data/mRNASeq/ISB025/ISB025_CT/Raw_data/FASTQ/' ,files$Folder[dir.row], '/', files$File[dir.row])
    cat(paste0(files.temp.out[row], '_R1_001.fastq.gz '))
  }
  
  out.sample <- files$Sample[dir.row]
  
  cat(paste0('> ', '/data/NCATS_ifx/data/mRNASeq/ISB025/ISB025_CT/Raw_data/concatenated_lanes/', out.sample, '/', out.sample, '_R1_001.fastq.gz'))
  cat("\n\n", sep="")
  
  rm(files.temp.out, files.row, files.temp, dir.row)
}
# R2----
for (sample in samples){
  files.row <- grep(sample, rows$Sample)
  files.temp <- unlist(as.numeric(str_split_fixed(rows$Rows[files.row], " ", Inf)))
  
  files.temp.out <- c()
  cat('cat ', sep='') 
  
  for (row in files.temp){
    files.temp[row] <- files$File[row]
    dir.row <- grep(files.temp[row],files$File)
    files.temp.out[row] <- paste0('/data/NCATS_ifx/data/mRNASeq/ISB025/ISB025_CT/Raw_data/FASTQ/' ,files$Folder[dir.row], '/', files$File[dir.row])
    cat(paste0(files.temp.out[row], '_R2_001.fastq.gz '))
  }
  
  out.sample <- files$Sample[dir.row]
  
  cat(paste0('> ', '/data/NCATS_ifx/data/mRNASeq/ISB025/ISB025_CT/Raw_data/concatenated_lanes/', out.sample, '/', out.sample, '_R2_001.fastq.gz'))
  cat("\n\n", sep="")
  
  rm(files.temp.out, files.row, files.temp, dir.row)
}
# swarm -f fastq_cat_081820.sh -t 1 -g 1 --time=8:00:00 #

# mapping-------
samples <- unique(files$Sample)

dir <- '/data/NCATS_ifx/data/mRNASeq/ISB025/ISB025_CT'
for (sample in samples){
  cat(paste0('cd ',dir,'/trimmomatic && java -jar $TRIMMOJAR PE -phred33 ', dir ,'/Raw_data/concatenated_lanes/',sample,'/', sample,'_R1_001.fastq.gz ', dir, '/Raw_data/concatenated_lanes/',sample, '/',
             sample, '_R2_001.fastq.gz ', '-baseout ', sample, '.fastq.gz ',
             'ILLUMINACLIP:/usr/local/apps/trimmomatic/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:10 TRAILING:5 MAXINFO:50:0.97 MINLEN:36 && cd ',dir,'/BAM && STAR --runThreadN $SLURM_CPUS_PER_TASK --genomeDir /fdb/STAR_current/GENCODE/Gencode_human/release_27/genes-150 --sjdbOverhang 150 --outSAMunmapped Within --outFilterType BySJout --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --sjdbScore 1 --readFilesIn ',dir,'/trimmomatic/',sample,'_1P.fastq.gz ',dir,'/trimmomatic/',sample,'_2P.fastq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ', sample,'_hg38 && htseq-count -f bam -r pos -s no -t exon -m union ',dir,'/BAM/',sample,'_hg38Aligned.sortedByCoord.out.bam /fdb/GENCODE/Gencode_human/release_27/gencode.v27.annotation.gtf --idattr gene_name > ',dir,'/htseq/',sample,'_htseq_counts.txt'), sep="\n")
  cat("\n")
}

# swarm -g 40 -t 12 --time=48:00:00 --module=trimmomatic,STAR,htseq -f ISB025_CT_pipeline.sh #