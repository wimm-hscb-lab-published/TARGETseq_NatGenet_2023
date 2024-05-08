# Introduction

This repository contains the scripts, metadata, and example output files for TARGETseq_2.0 as reported in Rodriguez-Meira _et al._ (Nature Genetics, 2023). Compared to its predecessor described in in Rodriguez-Meira _et al._ (Molecular Cell, 2019), the following bioinformatics-related improvements were implemented in TARGETseq_2.0:

1. Alignment of DNA-sequencing reads with Burrow-Wheeler Aligner (previously, DNA-sequencing reads were aligned with STAR aligner).

2. Detection of small-to-moderately sized insertions and deletions (indels) with Mutect2 software (previously, only had limited capacity for small indel detection with Samtools mpileup). Single-nucleotide variants (SNVs) remain detected with Samtools mpileup in TARGETseq_2.0.

3. Objective statistical approach for assigning genotypes into wildtype, heterozygous, or  homozygous mutants (previously, genotype assignments were made based solely on variant allele frequency(VAF)).

# Example data for exercise

Nine-hundred and fifty-four pairs (Read 1 and Read 2) of FASTQs are provided for an individual (anonymised ID: GST010) with one TP53 single-nucleotide variant (SNV) and one CALR 52-bp deletion. This provides us an opportunity to learn to process the data and call variants for the two main types of variants, i.e., SNV and indel.

Download link to FASTQs here: http://datashare.molbiol.ox.ac.uk/public/project/meadlab/wwen/Vladimir/GST010.tar.gz

When using your own set of FASTQs in the future, please ensure that one (pair) of FASTQ represents one cell, i.e., the original Illumina BCL files or FASTQ files should have been demultiplexed so that each (pair) of FASTQ file represents only one cell.

# Notes on scripts

- Script_02_DetectVariants_gDNA: Align sequencing reads, retrieve DNA-sequencing reads based on primer position, call variants, and assess coverage at expected variant sites. One file per cell.

- Script_02_DetectVariants_gDNA_perl: Need to be specified in Script_02_DetectVariants_gDNA. One file per cell.

- Script_03_DetectVariants_mRNA: Align sequencing reads, retrieve RNA-sequencing reads based on primer position, call variants, and assess coverage at expected variant sites. One file per cell.

- Script_03_DetectVariants_mRNA_perl: Need to be specified in Script_03_DetectVariants_mRNA. One file per cell.

- Script_04_1_MergeFiles_Coverage_gDNA.R: Collate coverage computed at variant sites.

- Script_04_2_MergeFiles_Mutect2_gDNA.R: Collate genotyped variants from Mutect2 output file.

- Script_04_3_MergeFiles_mpileup_gDNA.R: Collate genotyped variants from Samtools mpileup output file.

- Script_05_1-3: Same as Script_04_1-3 but for RNA-sequencing reads.

- Script_07_1_TabulateCounts_Mutect2_gDNA.R: Tabulate reference and mutant allele counts for genotyped indels as reported by Mutect2.

- Script_07_2_TabulateCounts_mpileup_gDNA.R: Tabulate reference and mutant allele counts for genotyped SNVs as reported by mpileup.

- Script_08_1-2: Same as Script_07_1-2 but for RNA-sequencing reads.

- Script_09_1_AnnotateGenotype_Mutect2_gDNA.R: Assign genotypes for indels.

- Script_09_2_AnnotateGenotype_mpileup_gDNA.R: Assign genotypes for SNVs.

- Script_10_1-2: Same as Script_09_1-2 but for RNA-sequencing reads.


# Notes on metadata

Located in the _metadata_ folder.

- SampleList_GST010.txt: Cell IDs for individual (anonymised ID: GST010). Helpful for retrieving this individual's FASTQs from the larger list of FASTQs provided.

- gPRIMERS_p53MPNAML.bed: List of primers used for genotyping DNA. This information will be used to retrieve DNA reads in Script_02_DetectVariants_gDNA_perl prior to downstream variant calling.

- mPRIMERS_p53MPNAML.bed: List of primers used for genotyping RNA. This information will be used to retrieve RNA reads in Script_03_DetectVariants_mRNA_perl prior to downstream variant calling.

- Intervals/bedtools/GST010_variantfile_NewGenoPipe.bed: Single-base position of genotyped variants. This information will be used by BEDTools to assess coverage at variant site. Variants with low coverage will be flagged in the final output files.

- Intervals/Mutect2/GST010_regionfile_NewGenoPipe.bed: Genomic range to call variants. Only applicable to indels. Typically specify plus/minus 200bp of the expected indel position. This tells Mutect2 to only look for variants specified in this range. If this file is not provided, Mutect2 will search the entire genome for variants, and this is computationally costly and inefficient.

- Intervals/mpileup/GST010_regionfile_NewGenoPipe.bed: Single genomic position to call variants. Only applicable to SNVs. This tells Samtools to only look for variants specified in at this specific genomic position.

- Variants_Metadata.xlsx: Genomic position, reference and mutant allele of genotyped variants.


# Notes on output file

Example output files for each step of the analysis in the _out_ folder. The ultimate files with the assigned genotypes are GST010_Mutect2_AlleleCounts_GenotypeAssigned_gDNA.txt and GST010_Mutect2_AlleleCounts_GenotypeAssigned_mRNA.txt

# Handy tips before getting started

- Always remember that SNVs are detected by Samtools mpileup whereas indels are detected by Mutect2. Although Mutect2 can in principle be used to detect SNVs, it often misses many (false negative).

- Note that in Variants_Metadata.xlsx, SNVs are typically straightforward to specify as they involved a single-base change. For indels, occasionally the reference and/or mutant alleles may be slightly different from the ones expected. If unsure, take a peak into the Mutect2 output file generated from Script_02_DetectVariants_gDNA and/or Script_03_DetectVariants_mRNA, and then proceed to populate the information in Variants_Metadata.xlsx. And then proceed with executing the rest of the scripts.

- Execute one script (step) at a time. We apologise for the lack of framework that integrates all steps into a single push-of-a-button.

- Ensure the relevant Linux and R packages are installed prior to executing the scripts.

- Ensure the paths are updated prior to executing the scripts.


# Contact

Adam Mead <adam.mead@imm.ox.ac.uk>

Sean Wen <sean.wen@astrazeneca.com> or <sw2136@cam.ac.uk>

Alba Rodriguez-Meira <alba_rodriguez-meira@dfci.harvard.edu> 

Affaf Aliouat <affaf.aliouat@ndcls.ox.ac.uk>
