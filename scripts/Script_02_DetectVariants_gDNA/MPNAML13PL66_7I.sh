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
mkdir -p /project/meadlab/wwen/Vladimir/GST010/MPNAML13PL66_7I/gDNA/
cd /project/meadlab/wwen/Vladimir/GST010/MPNAML13PL66_7I/gDNA/

###############################################################
##################### RETRIEVE gDNA READS #####################
###############################################################

# Align reads
bwa mem /project/meadlab/wwen/References/GRCh37_GENCODE_genome_sorted_bwa_indexed/GRCh37.primary_assembly.genome.fa \
        /t1-data/project/meadlab/arodrigu/MPNAML/scGenotyping/fastq/final/MPNAML13PL66_7I_R1.fq.gz \
        /t1-data/project/meadlab/arodrigu/MPNAML/scGenotyping/fastq/final/MPNAML13PL66_7I_R2.fq.gz > \
        MPNAML13PL66_7I_bwa_aligned.sam \
        -t 4

# Convert SAM to BAM
samtools view -bS -o MPNAML13PL66_7I_bwa_aligned.bam \
                     MPNAML13PL66_7I_bwa_aligned.sam \
                     -@ 4
	
# Sort by coordinates
samtools sort    MPNAML13PL66_7I_bwa_aligned.bam \
              -o MPNAML13PL66_7I.sorted.bam \
              -@ 4

# Index BAM
samtools index MPNAML13PL66_7I.sorted.bam

# Remove intermediate files
rm -rf MPNAML13PL66_7I_bwa_aligned.sam
rm -rf MPNAML13PL66_7I_bwa_aligned.bam

# Retrieve gDNA reads
perl /project/meadlab/wwen/Vladimir/Scripts/Scripts_01_GST010/Script_02_DetectVariants_gDNA_perl/MPNAML13PL66_7I.pl

###############################################################
#################### gDNA READS RETRIEVED #####################
###############################################################

# Concatenate BAM header to SAM header
samtools view -H MPNAML13PL66_7I.sorted.bam > MPNAML13PL66_7I_header.txt
          
cat MPNAML13PL66_7I_header.txt MPNAML13PL66_7I.sorted.genome.sam > MPNAML13PL66_7I.genome.withHeader.sam

# Convert SAM to BAM
samtools view -bS -o MPNAML13PL66_7I.genome.bam MPNAML13PL66_7I.genome.withHeader.sam -@ 4
       
# Sort by coordinates
samtools sort MPNAML13PL66_7I.genome.bam -o MPNAML13PL66_7I.sorted.bam -@ 4

# Index BAM
samtools index MPNAML13PL66_7I.sorted.bam

# Remove intermediate files
rm -rf MPNAML13PL66_7I_header.txt
rm -rf MPNAML13PL66_7I.sorted.genome.sam
rm -rf MPNAML13PL66_7I.genome.withHeader.sam
rm -rf MPNAML13PL66_7I.genome.bam

###############################################################
######################## CALL VARIANTS ########################
###############################################################

# Assign read groups
AddOrReplaceReadGroups CREATE_INDEX=True \
                       SORT_ORDER=coordinate \
                       RGPL=illumina \
                       RGSM=MPNAML13PL66_7I \
                       RGPU=MPNAML13PL66_7I \
                       RGLB=MPNAML13PL66_7I \
                       RGID=MPNAML13PL66_7I \
                       I=MPNAML13PL66_7I.sorted.bam \
                       O=MPNAML13PL66_7I.RGAssignedsorted.genome.bam

# Mark duplicates
MarkDuplicates CREATE_INDEX=True \
               I=MPNAML13PL66_7I.RGAssignedsorted.genome.bam \
               O=MPNAML13PL66_7I.RGAssignedsorted.DupMarked.genome.bam \
               M=MPNAML13PL66_7I.RGAssignedsorted.DupMarked.genome.metrics.txt
               
# Mutect2
gatk Mutect2 -R /project/meadlab/wwen/References/GRCh37_GENCODE_genome_sorted_bwa_indexed/GRCh37.primary_assembly.genome.fa \
             -L /project/meadlab/wwen/Vladimir/Metadata/Intervals/Mutect2/GST010_regionfile_NewGenoPipe.bed \
             -I MPNAML13PL66_7I.RGAssignedsorted.DupMarked.genome.bam \
             -O MPNAML13PL66_7I_MuTect2.vcf \
             --tumor-lod-to-emit 2.0 \
             --disable-read-filter NotDuplicateReadFilter \
             --max-reads-per-alignment-start 0 \
             --bam-output MPNAML13PL66_7I.haplotype.bam
        
# mpileup
samtools mpileup --count-orphans \
                 --max-depth 999999 \
                 --min-BQ 13 \
                 --ignore-overlaps \
                 --positions /project/meadlab/wwen/Vladimir/Metadata/Intervals/mpileup/GST010_regionfile_NewGenoPipe.bed \
                 --excl-flags SECONDARY \
                 --output MPNAML13PL66_7I_mpileup.txt \
                 MPNAML13PL66_7I.RGAssignedsorted.DupMarked.genome.bam

# Determine counts at variant site
bedtools coverage -counts \
                  -g /project/meadlab/wwen/References/GRCh37_GENCODE_genome_sorted_bwa_indexed/GRCh37.primary_assembly.genome.txt \
                  -sorted \
                  -a /project/meadlab/wwen/Vladimir/Metadata/Intervals/bedtools/GST010_variantfile_NewGenoPipe.bed \
                  -b MPNAML13PL66_7I.RGAssignedsorted.DupMarked.genome.bam > \
                     MPNAML13PL66_7I_Coverage.bed

# Remove intermediate files
rm -rf MPNAML13PL66_7I.sorted.bam
rm -rf MPNAML13PL66_7I.sorted.bam.bai
rm -rf MPNAML13PL66_7I.RGAssignedsorted.genome.bam
rm -rf MPNAML13PL66_7I.RGAssignedsorted.genome.bai
rm -rf MPNAML13PL66_7I.RGAssignedsorted.DupMarked.genome.metrics.txt
rm -rf MPNAML13PL66_7I_MuTect2.vcf.idx
rm -rf MPNAML13PL66_7I_MuTect2.vcf.stats
