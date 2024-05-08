#!/bin/sh

# Load modules
module load rna-star/2.6.0c
module load samtools/1.9
module load python-base/3.6.10
module load cbrg
module load picard-tools/2.3.0
module load gatk/4.2.0.0
module add bedtools/2.27.1

mkdir -p /project/meadlab/wwen/Vladimir/GST010/MPNAML13PL71_17B/mRNA/
cd /project/meadlab/wwen/Vladimir/GST010/MPNAML13PL71_17B/mRNA/

###############################################################
##################### RETRIEVE mRNA READS #####################
###############################################################

# Align reads
STAR --runThreadN 4 \
     --genomeDir  /project/meadlab/wwen/References/GRCh37_GENCODE_genome_STAR_indexed/ \
     --readFilesIn /t1-data/project/meadlab/arodrigu/MPNAML/scGenotyping/fastq/final/MPNAML13PL71_17B_R1.fq.gz /t1-data/project/meadlab/arodrigu/MPNAML/scGenotyping/fastq/final/MPNAML13PL71_17B_R2.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix /project/meadlab/wwen/Vladimir/GST010/MPNAML13PL71_17B/mRNA/MPNAML13PL71_17B. \
     --outTmpDir /project/meadlab/wwen/Vladimir/GST010/MPNAML13PL71_17B/mRNA/MPNAML13PL71_17B \
     --outReadsUnmapped Fastx \
     --limitBAMsortRAM 10000000000 \
     --outSAMtype BAM Unsorted

# Sort by coordinates
samtools sort MPNAML13PL71_17B.Aligned.out.bam \
          -o MPNAML13PL71_17B.sorted.bam \
          -@ 4

# Index BAM
samtools index MPNAML13PL71_17B.sorted.bam

# Remove intermediate files
rm -rf MPNAML13PL71_17B.Aligned.out.bam
rm -rf MPNAML13PL71_17B.Log.final.out
rm -rf MPNAML13PL71_17B.Log.out
rm -rf MPNAML13PL71_17B.Log.progress.out
rm -rf MPNAML13PL71_17B.SJ.out.tab
rm -rf MPNAML13PL71_17Btmp.fifo.read1
rm -rf MPNAML13PL71_17Btmp.fifo.read2
rm -rf MPNAML13PL71_17B.Unmapped.out.mate1
rm -rf MPNAML13PL71_17B.Unmapped.out.mate2

# Retrieve mRNA reads
perl /project/meadlab/wwen/Vladimir/Scripts/Scripts_01_GST010/Script_03_DetectVariants_mRNA_perl/MPNAML13PL71_17B.pl

###############################################################
#################### mRNA READS RETRIEVED #####################
###############################################################

# Concatenate BAM header to SAM header
samtools view -H MPNAML13PL71_17B.sorted.bam > MPNAML13PL71_17B_header.txt
          
cat MPNAML13PL71_17B_header.txt MPNAML13PL71_17B.sorted.genome.sam > MPNAML13PL71_17B.genome.withHeader.sam

# Convert SAM to BAM
samtools view -bS -o MPNAML13PL71_17B.genome.bam MPNAML13PL71_17B.genome.withHeader.sam -@ 4
       
# Sort by coordiates
samtools sort MPNAML13PL71_17B.genome.bam -o MPNAML13PL71_17B.sorted.bam -@ 4

# Index BAM
samtools index MPNAML13PL71_17B.sorted.bam

# Remove intermediate files
rm -rf MPNAML13PL71_17B_header.txt
rm -rf MPNAML13PL71_17B.sorted.genome.sam
rm -rf MPNAML13PL71_17B.genome.withHeader.sam
rm -rf MPNAML13PL71_17B.genome.bam

###############################################################
######################## CALL VARIANTS ########################
###############################################################

# Assign read groups
AddOrReplaceReadGroups CREATE_INDEX=True \
                   SORT_ORDER=coordinate \
                   RGPL=illumina \
                   RGSM=MPNAML13PL71_17B \
                   RGPU=MPNAML13PL71_17B \
                   RGLB=MPNAML13PL71_17B \
                   RGID=MPNAML13PL71_17B \
                   I=MPNAML13PL71_17B.sorted.bam \
                   O=MPNAML13PL71_17B.RGAssignedsorted.genome.bam

# Mark duplicates
MarkDuplicates CREATE_INDEX=True \
          I=MPNAML13PL71_17B.RGAssignedsorted.genome.bam \
          O=MPNAML13PL71_17B.RGAssignedsorted.DupMarked.genome.bam \
          M=MPNAML13PL71_17B.RGAssignedsorted.DupMarked.genome.metrics.txt

# Hard-clip splicing intervals
gatk SplitNCigarReads -R /project/meadlab/wwen/References/GRCh37_GENCODE_genome/GRCh37.primary_assembly.genome.fa \
                  -I MPNAML13PL71_17B.RGAssignedsorted.DupMarked.genome.bam \
                  -O MPNAML13PL71_17B.RGAssignedsorted.DupMark.HardClipped.bam

# MuTect2
gatk Mutect2 -R /project/meadlab/wwen/References/GRCh37_GENCODE_genome/GRCh37.primary_assembly.genome.fa \
             -L /project/meadlab/wwen/Vladimir/Metadata/Intervals/Mutect2/GST010_regionfile_NewGenoPipe.bed \
             -I MPNAML13PL71_17B.RGAssignedsorted.DupMark.HardClipped.bam \
             -O MPNAML13PL71_17B_MuTect2.vcf \
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
                 --output MPNAML13PL71_17B_mpileup.txt \
                 MPNAML13PL71_17B.RGAssignedsorted.DupMark.HardClipped.bam
                 
# Determine counts at variant site
bedtools coverage -counts \
             -g /project/meadlab/wwen/References/GRCh37_GENCODE_genome_sorted_bwa_indexed/GRCh37.primary_assembly.genome.txt \
             -sorted \
             -a /project/meadlab/wwen/Vladimir/Metadata/Intervals/bedtools/GST010_variantfile_NewGenoPipe.bed \
             -b MPNAML13PL71_17B.RGAssignedsorted.DupMark.HardClipped.bam > \
                MPNAML13PL71_17B_Coverage.bed

# Remove intermediate files
rm -rf MPNAML13PL71_17B_MuTect2.vcf.idx
rm -rf MPNAML13PL71_17B_MuTect2.vcf.stats
rm -rf MPNAML13PL71_17B.sorted.bam
rm -rf MPNAML13PL71_17B.sorted.bam.bai
rm -rf MPNAML13PL71_17B.RGAssignedsorted.genome.bam
rm -rf MPNAML13PL71_17B.RGAssignedsorted.genome.bai
rm -rf MPNAML13PL71_17B.RGAssignedsorted.DupMarked.genome.bam
rm -rf MPNAML13PL71_17B.RGAssignedsorted.DupMarked.genome.bai
rm -rf MPNAML13PL71_17B.RGAssignedsorted.DupMarked.genome.metrics.txt
