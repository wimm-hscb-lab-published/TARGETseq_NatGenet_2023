#!/bin/sh

# Load modules
module load bwa/0.7.17
module load samtools/1.9
module load python-base/3.6.10
module load cbrg
module load picard-tools/2.3.0
module load gatk/4.2.0.0
module add bedtools/2.27.1

# Create output directory
mkdir -p /project/meadlab/wwen/Vladimir/GST010/MPNAML13PL72_11P/gDNA/
cd /project/meadlab/wwen/Vladimir/GST010/MPNAML13PL72_11P/gDNA/

###############################################################
##################### RETRIEVE gDNA READS #####################
###############################################################

# Align reads
bwa mem /project/meadlab/wwen/References/GRCh37_GENCODE_genome_sorted_bwa_indexed/GRCh37.primary_assembly.genome.fa \
        /t1-data/project/meadlab/arodrigu/MPNAML/scGenotyping/fastq/final/MPNAML13PL72_11P_R1.fq.gz \
        /t1-data/project/meadlab/arodrigu/MPNAML/scGenotyping/fastq/final/MPNAML13PL72_11P_R2.fq.gz > \
        MPNAML13PL72_11P_bwa_aligned.sam \
        -t 4

# Convert SAM to BAM
samtools view -bS -o MPNAML13PL72_11P_bwa_aligned.bam \
                     MPNAML13PL72_11P_bwa_aligned.sam \
                     -@ 4
	
# Sort by coordinates
samtools sort    MPNAML13PL72_11P_bwa_aligned.bam \
              -o MPNAML13PL72_11P.sorted.bam \
              -@ 4

# Index BAM
samtools index MPNAML13PL72_11P.sorted.bam

# Remove intermediate files
rm -rf MPNAML13PL72_11P_bwa_aligned.sam
rm -rf MPNAML13PL72_11P_bwa_aligned.bam

# Retrieve gDNA reads
perl /project/meadlab/wwen/Vladimir/Scripts/Scripts_01_GST010/Script_02_DetectVariants_gDNA_perl/MPNAML13PL72_11P.pl

###############################################################
#################### gDNA READS RETRIEVED #####################
###############################################################

# Concatenate BAM header to SAM header
samtools view -H MPNAML13PL72_11P.sorted.bam > MPNAML13PL72_11P_header.txt
          
cat MPNAML13PL72_11P_header.txt MPNAML13PL72_11P.sorted.genome.sam > MPNAML13PL72_11P.genome.withHeader.sam

# Convert SAM to BAM
samtools view -bS -o MPNAML13PL72_11P.genome.bam MPNAML13PL72_11P.genome.withHeader.sam -@ 4
       
# Sort by coordinates
samtools sort MPNAML13PL72_11P.genome.bam -o MPNAML13PL72_11P.sorted.bam -@ 4

# Index BAM
samtools index MPNAML13PL72_11P.sorted.bam

# Remove intermediate files
rm -rf MPNAML13PL72_11P_header.txt
rm -rf MPNAML13PL72_11P.sorted.genome.sam
rm -rf MPNAML13PL72_11P.genome.withHeader.sam
rm -rf MPNAML13PL72_11P.genome.bam

###############################################################
######################## CALL VARIANTS ########################
###############################################################

# Assign read groups
AddOrReplaceReadGroups CREATE_INDEX=True \
                       SORT_ORDER=coordinate \
                       RGPL=illumina \
                       RGSM=MPNAML13PL72_11P \
                       RGPU=MPNAML13PL72_11P \
                       RGLB=MPNAML13PL72_11P \
                       RGID=MPNAML13PL72_11P \
                       I=MPNAML13PL72_11P.sorted.bam \
                       O=MPNAML13PL72_11P.RGAssignedsorted.genome.bam

# Mark duplicates
MarkDuplicates CREATE_INDEX=True \
               I=MPNAML13PL72_11P.RGAssignedsorted.genome.bam \
               O=MPNAML13PL72_11P.RGAssignedsorted.DupMarked.genome.bam \
               M=MPNAML13PL72_11P.RGAssignedsorted.DupMarked.genome.metrics.txt
               
# Mutect2
gatk Mutect2 -R /project/meadlab/wwen/References/GRCh37_GENCODE_genome_sorted_bwa_indexed/GRCh37.primary_assembly.genome.fa \
             -L /project/meadlab/wwen/Vladimir/Metadata/Intervals/Mutect2/GST010_regionfile_NewGenoPipe.bed \
             -I MPNAML13PL72_11P.RGAssignedsorted.DupMarked.genome.bam \
             -O MPNAML13PL72_11P_MuTect2.vcf \
             --tumor-lod-to-emit 2.0 \
             --disable-read-filter NotDuplicateReadFilter \
             --max-reads-per-alignment-start 0 \
             --bam-output MPNAML13PL72_11P.haplotype.bam
        
# mpileup
samtools mpileup --count-orphans \
                 --max-depth 999999 \
                 --min-BQ 13 \
                 --ignore-overlaps \
                 --positions /project/meadlab/wwen/Vladimir/Metadata/Intervals/mpileup/GST010_regionfile_NewGenoPipe.bed \
                 --excl-flags SECONDARY \
                 --output MPNAML13PL72_11P_mpileup.txt \
                 MPNAML13PL72_11P.RGAssignedsorted.DupMarked.genome.bam

# Determine counts at variant site
bedtools coverage -counts \
                  -g /project/meadlab/wwen/References/GRCh37_GENCODE_genome_sorted_bwa_indexed/GRCh37.primary_assembly.genome.txt \
                  -sorted \
                  -a /project/meadlab/wwen/Vladimir/Metadata/Intervals/bedtools/GST010_variantfile_NewGenoPipe.bed \
                  -b MPNAML13PL72_11P.RGAssignedsorted.DupMarked.genome.bam > \
                     MPNAML13PL72_11P_Coverage.bed

# Remove intermediate files
rm -rf MPNAML13PL72_11P.sorted.bam
rm -rf MPNAML13PL72_11P.sorted.bam.bai
rm -rf MPNAML13PL72_11P.RGAssignedsorted.genome.bam
rm -rf MPNAML13PL72_11P.RGAssignedsorted.genome.bai
rm -rf MPNAML13PL72_11P.RGAssignedsorted.DupMarked.genome.metrics.txt
rm -rf MPNAML13PL72_11P_MuTect2.vcf.idx
rm -rf MPNAML13PL72_11P_MuTect2.vcf.stats
