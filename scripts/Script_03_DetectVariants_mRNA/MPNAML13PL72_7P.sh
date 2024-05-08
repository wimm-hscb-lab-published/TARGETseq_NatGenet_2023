#!/bin/sh

# Load modules
module load rna-star/2.6.0c
module load samtools/1.9
module load python-base/3.6.10
module load cbrg
module load picard-tools/2.3.0
module load gatk/4.2.0.0
module add bedtools/2.27.1

mkdir -p /project/meadlab/wwen/Vladimir/GST010/MPNAML13PL72_7P/mRNA/
cd /project/meadlab/wwen/Vladimir/GST010/MPNAML13PL72_7P/mRNA/

###############################################################
##################### RETRIEVE mRNA READS #####################
###############################################################

# Align reads
STAR --runThreadN 4 \
     --genomeDir  /project/meadlab/wwen/References/GRCh37_GENCODE_genome_STAR_indexed/ \
     --readFilesIn /t1-data/project/meadlab/arodrigu/MPNAML/scGenotyping/fastq/final/MPNAML13PL72_7P_R1.fq.gz /t1-data/project/meadlab/arodrigu/MPNAML/scGenotyping/fastq/final/MPNAML13PL72_7P_R2.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix /project/meadlab/wwen/Vladimir/GST010/MPNAML13PL72_7P/mRNA/MPNAML13PL72_7P. \
     --outTmpDir /project/meadlab/wwen/Vladimir/GST010/MPNAML13PL72_7P/mRNA/MPNAML13PL72_7P \
     --outReadsUnmapped Fastx \
     --limitBAMsortRAM 10000000000 \
     --outSAMtype BAM Unsorted

# Sort by coordinates
samtools sort MPNAML13PL72_7P.Aligned.out.bam \
          -o MPNAML13PL72_7P.sorted.bam \
          -@ 4

# Index BAM
samtools index MPNAML13PL72_7P.sorted.bam

# Remove intermediate files
rm -rf MPNAML13PL72_7P.Aligned.out.bam
rm -rf MPNAML13PL72_7P.Log.final.out
rm -rf MPNAML13PL72_7P.Log.out
rm -rf MPNAML13PL72_7P.Log.progress.out
rm -rf MPNAML13PL72_7P.SJ.out.tab
rm -rf MPNAML13PL72_7Ptmp.fifo.read1
rm -rf MPNAML13PL72_7Ptmp.fifo.read2
rm -rf MPNAML13PL72_7P.Unmapped.out.mate1
rm -rf MPNAML13PL72_7P.Unmapped.out.mate2

# Retrieve mRNA reads
perl /project/meadlab/wwen/Vladimir/Scripts/Scripts_01_GST010/Script_03_DetectVariants_mRNA_perl/MPNAML13PL72_7P.pl

###############################################################
#################### mRNA READS RETRIEVED #####################
###############################################################

# Concatenate BAM header to SAM header
samtools view -H MPNAML13PL72_7P.sorted.bam > MPNAML13PL72_7P_header.txt
          
cat MPNAML13PL72_7P_header.txt MPNAML13PL72_7P.sorted.genome.sam > MPNAML13PL72_7P.genome.withHeader.sam

# Convert SAM to BAM
samtools view -bS -o MPNAML13PL72_7P.genome.bam MPNAML13PL72_7P.genome.withHeader.sam -@ 4
       
# Sort by coordiates
samtools sort MPNAML13PL72_7P.genome.bam -o MPNAML13PL72_7P.sorted.bam -@ 4

# Index BAM
samtools index MPNAML13PL72_7P.sorted.bam

# Remove intermediate files
rm -rf MPNAML13PL72_7P_header.txt
rm -rf MPNAML13PL72_7P.sorted.genome.sam
rm -rf MPNAML13PL72_7P.genome.withHeader.sam
rm -rf MPNAML13PL72_7P.genome.bam

###############################################################
######################## CALL VARIANTS ########################
###############################################################

# Assign read groups
AddOrReplaceReadGroups CREATE_INDEX=True \
                   SORT_ORDER=coordinate \
                   RGPL=illumina \
                   RGSM=MPNAML13PL72_7P \
                   RGPU=MPNAML13PL72_7P \
                   RGLB=MPNAML13PL72_7P \
                   RGID=MPNAML13PL72_7P \
                   I=MPNAML13PL72_7P.sorted.bam \
                   O=MPNAML13PL72_7P.RGAssignedsorted.genome.bam

# Mark duplicates
MarkDuplicates CREATE_INDEX=True \
          I=MPNAML13PL72_7P.RGAssignedsorted.genome.bam \
          O=MPNAML13PL72_7P.RGAssignedsorted.DupMarked.genome.bam \
          M=MPNAML13PL72_7P.RGAssignedsorted.DupMarked.genome.metrics.txt

# Hard-clip splicing intervals
gatk SplitNCigarReads -R /project/meadlab/wwen/References/GRCh37_GENCODE_genome/GRCh37.primary_assembly.genome.fa \
                  -I MPNAML13PL72_7P.RGAssignedsorted.DupMarked.genome.bam \
                  -O MPNAML13PL72_7P.RGAssignedsorted.DupMark.HardClipped.bam

# MuTect2
gatk Mutect2 -R /project/meadlab/wwen/References/GRCh37_GENCODE_genome/GRCh37.primary_assembly.genome.fa \
             -L /project/meadlab/wwen/Vladimir/Metadata/Intervals/Mutect2/GST010_regionfile_NewGenoPipe.bed \
             -I MPNAML13PL72_7P.RGAssignedsorted.DupMark.HardClipped.bam \
             -O MPNAML13PL72_7P_MuTect2.vcf \
             --tumor-lod-to-emit 2.0 \
             --disable-read-filter NotDuplicateReadFilter \
             --max-reads-per-alignment-start 0

# mpileup
samtools mpileup --count-orphans \
                 --max-depth 999999 \
                 --min-BQ 13 \
                 --ignore-overlaps \
                 --positions /project/meadlab/wwen/Vladimir/Metadata/Intervals/mpileup/GST010_regionfile_NewGenoPipe.bed \
                 --excl-flags SECONDARY \
                 --output MPNAML13PL72_7P_mpileup.txt \
                 MPNAML13PL72_7P.RGAssignedsorted.DupMark.HardClipped.bam
                 
# Determine counts at variant site
bedtools coverage -counts \
             -g /project/meadlab/wwen/References/GRCh37_GENCODE_genome_sorted_bwa_indexed/GRCh37.primary_assembly.genome.txt \
             -sorted \
             -a /project/meadlab/wwen/Vladimir/Metadata/Intervals/bedtools/GST010_variantfile_NewGenoPipe.bed \
             -b MPNAML13PL72_7P.RGAssignedsorted.DupMark.HardClipped.bam > \
                MPNAML13PL72_7P_Coverage.bed

# Remove intermediate files
rm -rf MPNAML13PL72_7P_MuTect2.vcf.idx
rm -rf MPNAML13PL72_7P_MuTect2.vcf.stats
rm -rf MPNAML13PL72_7P.sorted.bam
rm -rf MPNAML13PL72_7P.sorted.bam.bai
rm -rf MPNAML13PL72_7P.RGAssignedsorted.genome.bam
rm -rf MPNAML13PL72_7P.RGAssignedsorted.genome.bai
rm -rf MPNAML13PL72_7P.RGAssignedsorted.DupMarked.genome.bam
rm -rf MPNAML13PL72_7P.RGAssignedsorted.DupMarked.genome.bai
rm -rf MPNAML13PL72_7P.RGAssignedsorted.DupMarked.genome.metrics.txt