#!/bin/sh

# Load modules
module load rna-star/2.6.0c
module load samtools/1.9
module load python-base/3.6.10
module load cbrg
module load picard-tools/2.3.0
module load gatk/4.2.0.0
module add bedtools/2.27.1

mkdir -p /project/meadlab/wwen/Vladimir/GST010/MPNAML13PL64_23O/mRNA/
cd /project/meadlab/wwen/Vladimir/GST010/MPNAML13PL64_23O/mRNA/

###############################################################
##################### RETRIEVE mRNA READS #####################
###############################################################

# Align reads
STAR --runThreadN 4 \
     --genomeDir  /project/meadlab/wwen/References/GRCh37_GENCODE_genome_STAR_indexed/ \
     --readFilesIn /t1-data/project/meadlab/arodrigu/MPNAML/scGenotyping/fastq/final/MPNAML13PL64_23O_R1.fq.gz /t1-data/project/meadlab/arodrigu/MPNAML/scGenotyping/fastq/final/MPNAML13PL64_23O_R2.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix /project/meadlab/wwen/Vladimir/GST010/MPNAML13PL64_23O/mRNA/MPNAML13PL64_23O. \
     --outTmpDir /project/meadlab/wwen/Vladimir/GST010/MPNAML13PL64_23O/mRNA/MPNAML13PL64_23O \
     --outReadsUnmapped Fastx \
     --limitBAMsortRAM 10000000000 \
     --outSAMtype BAM Unsorted

# Sort by coordinates
samtools sort MPNAML13PL64_23O.Aligned.out.bam \
          -o MPNAML13PL64_23O.sorted.bam \
          -@ 4

# Index BAM
samtools index MPNAML13PL64_23O.sorted.bam

# Remove intermediate files
rm -rf MPNAML13PL64_23O.Aligned.out.bam
rm -rf MPNAML13PL64_23O.Log.final.out
rm -rf MPNAML13PL64_23O.Log.out
rm -rf MPNAML13PL64_23O.Log.progress.out
rm -rf MPNAML13PL64_23O.SJ.out.tab
rm -rf MPNAML13PL64_23Otmp.fifo.read1
rm -rf MPNAML13PL64_23Otmp.fifo.read2
rm -rf MPNAML13PL64_23O.Unmapped.out.mate1
rm -rf MPNAML13PL64_23O.Unmapped.out.mate2

# Retrieve mRNA reads
perl /project/meadlab/wwen/Vladimir/Scripts/Scripts_01_GST010/Script_03_DetectVariants_mRNA_perl/MPNAML13PL64_23O.pl

###############################################################
#################### mRNA READS RETRIEVED #####################
###############################################################

# Concatenate BAM header to SAM header
samtools view -H MPNAML13PL64_23O.sorted.bam > MPNAML13PL64_23O_header.txt
          
cat MPNAML13PL64_23O_header.txt MPNAML13PL64_23O.sorted.genome.sam > MPNAML13PL64_23O.genome.withHeader.sam

# Convert SAM to BAM
samtools view -bS -o MPNAML13PL64_23O.genome.bam MPNAML13PL64_23O.genome.withHeader.sam -@ 4
       
# Sort by coordiates
samtools sort MPNAML13PL64_23O.genome.bam -o MPNAML13PL64_23O.sorted.bam -@ 4

# Index BAM
samtools index MPNAML13PL64_23O.sorted.bam

# Remove intermediate files
rm -rf MPNAML13PL64_23O_header.txt
rm -rf MPNAML13PL64_23O.sorted.genome.sam
rm -rf MPNAML13PL64_23O.genome.withHeader.sam
rm -rf MPNAML13PL64_23O.genome.bam

###############################################################
######################## CALL VARIANTS ########################
###############################################################

# Assign read groups
AddOrReplaceReadGroups CREATE_INDEX=True \
                   SORT_ORDER=coordinate \
                   RGPL=illumina \
                   RGSM=MPNAML13PL64_23O \
                   RGPU=MPNAML13PL64_23O \
                   RGLB=MPNAML13PL64_23O \
                   RGID=MPNAML13PL64_23O \
                   I=MPNAML13PL64_23O.sorted.bam \
                   O=MPNAML13PL64_23O.RGAssignedsorted.genome.bam

# Mark duplicates
MarkDuplicates CREATE_INDEX=True \
          I=MPNAML13PL64_23O.RGAssignedsorted.genome.bam \
          O=MPNAML13PL64_23O.RGAssignedsorted.DupMarked.genome.bam \
          M=MPNAML13PL64_23O.RGAssignedsorted.DupMarked.genome.metrics.txt

# Hard-clip splicing intervals
gatk SplitNCigarReads -R /project/meadlab/wwen/References/GRCh37_GENCODE_genome/GRCh37.primary_assembly.genome.fa \
                  -I MPNAML13PL64_23O.RGAssignedsorted.DupMarked.genome.bam \
                  -O MPNAML13PL64_23O.RGAssignedsorted.DupMark.HardClipped.bam

# MuTect2
gatk Mutect2 -R /project/meadlab/wwen/References/GRCh37_GENCODE_genome/GRCh37.primary_assembly.genome.fa \
             -L /project/meadlab/wwen/Vladimir/Metadata/Intervals/Mutect2/GST010_regionfile_NewGenoPipe.bed \
             -I MPNAML13PL64_23O.RGAssignedsorted.DupMark.HardClipped.bam \
             -O MPNAML13PL64_23O_MuTect2.vcf \
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
                 --output MPNAML13PL64_23O_mpileup.txt \
                 MPNAML13PL64_23O.RGAssignedsorted.DupMark.HardClipped.bam
                 
# Determine counts at variant site
bedtools coverage -counts \
             -g /project/meadlab/wwen/References/GRCh37_GENCODE_genome_sorted_bwa_indexed/GRCh37.primary_assembly.genome.txt \
             -sorted \
             -a /project/meadlab/wwen/Vladimir/Metadata/Intervals/bedtools/GST010_variantfile_NewGenoPipe.bed \
             -b MPNAML13PL64_23O.RGAssignedsorted.DupMark.HardClipped.bam > \
                MPNAML13PL64_23O_Coverage.bed

# Remove intermediate files
rm -rf MPNAML13PL64_23O_MuTect2.vcf.idx
rm -rf MPNAML13PL64_23O_MuTect2.vcf.stats
rm -rf MPNAML13PL64_23O.sorted.bam
rm -rf MPNAML13PL64_23O.sorted.bam.bai
rm -rf MPNAML13PL64_23O.RGAssignedsorted.genome.bam
rm -rf MPNAML13PL64_23O.RGAssignedsorted.genome.bai
rm -rf MPNAML13PL64_23O.RGAssignedsorted.DupMarked.genome.bam
rm -rf MPNAML13PL64_23O.RGAssignedsorted.DupMarked.genome.bai
rm -rf MPNAML13PL64_23O.RGAssignedsorted.DupMarked.genome.metrics.txt