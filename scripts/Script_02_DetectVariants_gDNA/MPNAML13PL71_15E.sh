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
mkdir -p /project/meadlab/wwen/Vladimir/GST010/MPNAML13PL71_15E/gDNA/
cd /project/meadlab/wwen/Vladimir/GST010/MPNAML13PL71_15E/gDNA/

###############################################################
##################### RETRIEVE gDNA READS #####################
###############################################################

# Align reads
bwa mem /project/meadlab/wwen/References/GRCh37_GENCODE_genome_sorted_bwa_indexed/GRCh37.primary_assembly.genome.fa \
        /t1-data/project/meadlab/arodrigu/MPNAML/scGenotyping/fastq/final/MPNAML13PL71_15E_R1.fq.gz \
        /t1-data/project/meadlab/arodrigu/MPNAML/scGenotyping/fastq/final/MPNAML13PL71_15E_R2.fq.gz > \
        MPNAML13PL71_15E_bwa_aligned.sam \
        -t 4

# Convert SAM to BAM
samtools view -bS -o MPNAML13PL71_15E_bwa_aligned.bam \
                     MPNAML13PL71_15E_bwa_aligned.sam \
                     -@ 4
	
# Sort by coordinates
samtools sort    MPNAML13PL71_15E_bwa_aligned.bam \
              -o MPNAML13PL71_15E.sorted.bam \
              -@ 4

# Index BAM
samtools index MPNAML13PL71_15E.sorted.bam

# Remove intermediate files
rm -rf MPNAML13PL71_15E_bwa_aligned.sam
rm -rf MPNAML13PL71_15E_bwa_aligned.bam

# Retrieve gDNA reads
perl /project/meadlab/wwen/Vladimir/Scripts/Scripts_01_GST010/Script_02_DetectVariants_gDNA_perl/MPNAML13PL71_15E.pl

###############################################################
#################### gDNA READS RETRIEVED #####################
###############################################################

# Concatenate BAM header to SAM header
samtools view -H MPNAML13PL71_15E.sorted.bam > MPNAML13PL71_15E_header.txt
          
cat MPNAML13PL71_15E_header.txt MPNAML13PL71_15E.sorted.genome.sam > MPNAML13PL71_15E.genome.withHeader.sam

# Convert SAM to BAM
samtools view -bS -o MPNAML13PL71_15E.genome.bam MPNAML13PL71_15E.genome.withHeader.sam -@ 4
       
# Sort by coordinates
samtools sort MPNAML13PL71_15E.genome.bam -o MPNAML13PL71_15E.sorted.bam -@ 4

# Index BAM
samtools index MPNAML13PL71_15E.sorted.bam

# Remove intermediate files
rm -rf MPNAML13PL71_15E_header.txt
rm -rf MPNAML13PL71_15E.sorted.genome.sam
rm -rf MPNAML13PL71_15E.genome.withHeader.sam
rm -rf MPNAML13PL71_15E.genome.bam

###############################################################
######################## CALL VARIANTS ########################
###############################################################

# Assign read groups
AddOrReplaceReadGroups CREATE_INDEX=True \
                       SORT_ORDER=coordinate \
                       RGPL=illumina \
                       RGSM=MPNAML13PL71_15E \
                       RGPU=MPNAML13PL71_15E \
                       RGLB=MPNAML13PL71_15E \
                       RGID=MPNAML13PL71_15E \
                       I=MPNAML13PL71_15E.sorted.bam \
                       O=MPNAML13PL71_15E.RGAssignedsorted.genome.bam

# Mark duplicates
MarkDuplicates CREATE_INDEX=True \
               I=MPNAML13PL71_15E.RGAssignedsorted.genome.bam \
               O=MPNAML13PL71_15E.RGAssignedsorted.DupMarked.genome.bam \
               M=MPNAML13PL71_15E.RGAssignedsorted.DupMarked.genome.metrics.txt
               
# Mutect2
gatk Mutect2 -R /project/meadlab/wwen/References/GRCh37_GENCODE_genome_sorted_bwa_indexed/GRCh37.primary_assembly.genome.fa \
             -L /project/meadlab/wwen/Vladimir/Metadata/Intervals/Mutect2/GST010_regionfile_NewGenoPipe.bed \
             -I MPNAML13PL71_15E.RGAssignedsorted.DupMarked.genome.bam \
             -O MPNAML13PL71_15E_MuTect2.vcf \
             --tumor-lod-to-emit 2.0 \
             --disable-read-filter NotDuplicateReadFilter \
             --max-reads-per-alignment-start 0 \
             --bam-output MPNAML13PL71_15E.haplotype.bam
        
# mpileup
samtools mpileup --count-orphans \
                 --max-depth 999999 \
                 --min-BQ 13 \
                 --ignore-overlaps \
                 --positions /project/meadlab/wwen/Vladimir/Metadata/Intervals/mpileup/GST010_regionfile_NewGenoPipe.bed \
                 --excl-flags SECONDARY \
                 --output MPNAML13PL71_15E_mpileup.txt \
                 MPNAML13PL71_15E.RGAssignedsorted.DupMarked.genome.bam

# Determine counts at variant site
bedtools coverage -counts \
                  -g /project/meadlab/wwen/References/GRCh37_GENCODE_genome_sorted_bwa_indexed/GRCh37.primary_assembly.genome.txt \
                  -sorted \
                  -a /project/meadlab/wwen/Vladimir/Metadata/Intervals/bedtools/GST010_variantfile_NewGenoPipe.bed \
                  -b MPNAML13PL71_15E.RGAssignedsorted.DupMarked.genome.bam > \
                     MPNAML13PL71_15E_Coverage.bed

# Remove intermediate files
rm -rf MPNAML13PL71_15E.sorted.bam
rm -rf MPNAML13PL71_15E.sorted.bam.bai
rm -rf MPNAML13PL71_15E.RGAssignedsorted.genome.bam
rm -rf MPNAML13PL71_15E.RGAssignedsorted.genome.bai
rm -rf MPNAML13PL71_15E.RGAssignedsorted.DupMarked.genome.metrics.txt
rm -rf MPNAML13PL71_15E_MuTect2.vcf.idx
rm -rf MPNAML13PL71_15E_MuTect2.vcf.stats