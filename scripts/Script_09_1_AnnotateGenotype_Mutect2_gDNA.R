# Load packages
library(data.table)
library(gdata)
library(plyr)
library(reshape2)

# Read files
    # Allele counts
    path <- "/Users/seanwen/Documents/Vladimir/out/"
    file <- "GST010_Mutect2_AlleleCounts_gDNA.txt"
    df <- as.data.frame(fread(paste(path, file, sep=""), sep="\t", header=TRUE, stringsAsFactors=FALSE))
    
    # Variant metadata
    path <- "/Users/seanwen/Documents/Vladimir/Metadata/"
    file <- "Variants_Metadata.xlsx"
    md.variant <- read.xls(paste(path, file, sep=""), sheet=1, header=TRUE, stringsAsFactors=FALSE)
    
    # Coverage
    path <- "/Users/seanwen/Documents/Vladimir/out/"
    file <- "GST010_Coverage_gDNA.txt"
    cov <- as.data.frame(fread(paste(path, file, sep=""), sep="\t", header=TRUE, stringsAsFactors=FALSE))

# Subset INDELs only
variant_types <- c("Deletion", "Insertion")
md.variant <- md.variant[which(md.variant$variant_type %in% variant_types), ]

############################################################
##################### ANNOTATE COVERAGE ####################
############################################################

# Define variants
variant_names <- unique(df$variant_name)

.list <- list()

for(i in 1:length(variant_names)) {

    # Subset variant
    cov.small <- cov[which(cov$variant_name == variant_names[i]), ]
    df.small <- df[which(df$variant_name == variant_names[i]), ]
    md.variant.small <- md.variant[which(md.variant$variant_name == variant_names[i]), ]
    
    # Annotate coverage
    results <- join(df.small, cov.small[,c("cell.id", "coverage", "variant_name")], by="cell.id", type="full")
    
    # Annotate variant type for cells with no VCF records
    results$variant_type[is.na(results$variant_type)] <- md.variant.small$variant_type
    
    
    # Save into list
    .list[[i]] <- results

}

df <- do.call(rbind.data.frame, .list)

############################################################
################### ASSIGN GENOTYPE: INS ###################
############################################################

# Test
#names(df)[which(names(df)=="ins")] <- "ins.temp"
#names(df)[which(names(df)=="del")] <- "del.temp"
#names(df)[which(names(df)=="ins.temp")] <- "del"
#names(df)[which(names(df)=="del.temp")] <- "ins"
#df$variant_type[which(df$variant_type=="Deletion")] <- "Insertion"

# Subset variant type
table(df$variant_type)
df.small <- df[which(df$variant_type=="Insertion"), ]

.list.ins <- list()

# Assign genotype: VCF records + sufficient coverage
df.small. <- df.small[which(!is.na(df.small$chr)), ]
#df.small. <- df.small.[which(df.small.$coverage > df.small.$coverage.blank), ]
df.small. <- df.small.[which(df.small.$coverage > 10), ]

if(nrow(df.small.) != 0) {

    variant_names <- unique(df.small.$variant_name)

    .list <- list()

    for(k in 1:length(variant_names)) {

        # Subset variant
        df.small.. <- df.small.[which(df.small.$variant_name==variant_names[k]), ]
        
        if(nrow(df.small..) != 0) {
        
            # Retrieve ref, alt nt
            nt.ref <- unique(md.variant[which(md.variant$variant_name==variant_names[k]), "ref"])
            nt.ref <- "A"
            nt.alt <- "ins"
                
            # Retrieve ref, alt counts
            count.ref <- df.small..[, nt.ref]
            count.alt <- df.small..[, nt.alt]
            
            # Compute VAF
            df.small..$vaf <- count.alt / (count.ref + count.alt) * 100

            # Compute genotype score for each sample
            scores.list <- list()

            for(j in 1:length(count.ref)) {

                # Retrieve counts
                counts.wt.cell <- count.ref[j]
                counts.mt.cell <- count.alt[j]
                
                if(counts.wt.cell==0 & counts.mt.cell==0) {
                
                    counts.wt.cell <- 1
                    counts.mt.cell <- 0
                
                }
                
                # Define probabilities for each model
                prob.wt <- c(0.999, 0.5, 0.001)
                prob.mt <- c(0.001, 0.5, 0.999)
                
                scores <- NULL
                
                # Test each model
                for(i in 1:length(prob.wt)) {
                
                    # Retrieve test statistic
                    x2 <- chisq.test(c(counts.wt.cell, counts.mt.cell), p=c(prob.wt[i], prob.mt[i]))$statistic
                    
                    # Convert to zygosity score and save into vector
                    scores[i] <- 1/log10(x2 + 1)
                    
                }
                
                # Save into list
                scores.list[[j]] <- scores
                
                #print(j)
                    
            }

            # Tabulate results
            results <- do.call(rbind.data.frame, scores.list)

            # Assign genotype
            . <- apply(results, 1, function(x) {which.max(x)})
            .[which(.== "1")] <- "WT"
            .[which(.== "2")] <- "Het"
            .[which(.== "3")] <- "Hom"
            df.small..$genotype <- .
            
            # Save into list
            .list[[k]] <- df.small..
            
        } else {
        
            # Save into list
            .list[[k]] <- NULL
        
        }
        
    }

    results <- do.call(rbind.data.frame, .list)
    print(nrow(results)==nrow(df.small.))
    .list.ins[[1]] <- results
    
}
    
# Assign genotype: VCF records + low coverage
df.small. <- df.small[which(!is.na(df.small$chr)), ]
#df.small. <- df.small.[which(df.small.$coverage <= df.small.$coverage.blank), ]
df.small. <- df.small.[which(df.small.$coverage <= 10), ]

if(nrow(df.small.) != 0) {

    .list <- list()

    variant_names <- unique(df.small.$variant_name)
    
    for(k in 1:length(variant_names)) {

        # Subset variant
        df.small.. <- df.small.[which(df.small.$variant_name==variant_names[k]), ]
        
        if(nrow(df.small..) != 0) {
                
            # Retrieve ref, alt nt
            nt.ref <- unique(md.variant[which(md.variant$variant_name==variant_names[k]), "ref"])
            nt.ref <- "A"
            nt.alt <- "ins"
                
            # Retrieve ref, alt counts
            count.ref <- df.small..[, nt.ref]
            count.alt <- df.small..[, nt.alt]
            
            # Compute VAF
            df.small..$vaf <- count.alt / (count.ref + count.alt) * 100
            
            # Save into list
            .list[[k]] <- df.small..
            
        } else {
        
            # Save into list
            .list[[k]] <- NULL
        
        }

    }

    results <- do.call(rbind.data.frame, .list)
    print(nrow(results)==nrow(df.small.))
    df.small. <- results
    df.small.$genotype <- "Low Coverage"
    .list.ins[[2]] <- df.small.
    
}
    
# Assign genotype: No VCF records
df.small. <- df.small[which(is.na(df.small$chr)), ]

if(nrow(df.small.) != 0) {
    
    df.small.$vaf <- 0
    #df.small.$genotype <- ifelse(df.small.$coverage > df.small.$coverage.blank, "WT", "Low Coverage")
    df.small.$genotype <- ifelse(df.small.$coverage > 10, "WT", "Low Coverage")
    .list.ins[[3]] <- df.small.

}

# Merge
#sum(sapply(.list.ins, function(x) {nrow(x)})) == nrow(df.small)
df.ins <- do.call(rbind.data.frame, .list.ins)

############################################################
################### ASSIGN GENOTYPE: DEL ###################
############################################################

# Subset variant type
table(df$variant_type)
df.small <- df[which(df$variant_type=="Deletion"), ]

.list.del <- list()

# Assign genotype: VCF records + sufficient coverage
df.small. <- df.small[which(!is.na(df.small$chr)), ]
#df.small. <- df.small.[which(df.small.$coverage > df.small.$coverage.blank), ]
df.small. <- df.small.[which(df.small.$coverage > 10), ]

if(nrow(df.small.) != 0) {

    variant_names <- unique(df.small.$variant_name)

    .list <- list()

    for(k in 1:length(variant_names)) {

        # Subset variant
        df.small.. <- df.small.[which(df.small.$variant_name==variant_names[k]), ]
        
        if(nrow(df.small..) != 0) {
                
            # Retrieve ref, alt nt
            nt.ref <- unique(md.variant[which(md.variant$variant_name==variant_names[k]), "alt"])
            nt.alt <- "del"
                
            # Retrieve ref, alt counts
            count.ref <- df.small..[, nt.ref]
            count.alt <- df.small..[, nt.alt]
            
            # Compute VAF
            df.small..$vaf <- count.alt / (count.ref + count.alt) * 100

            # Compute genotype score for each sample
            scores.list <- list()

            for(j in 1:length(count.ref)) {

                # Retrieve counts
                counts.wt.cell <- count.ref[j]
                counts.mt.cell <- count.alt[j]
                
                if(counts.wt.cell==0 & counts.mt.cell==0) {
                
                    counts.wt.cell <- 1
                    counts.mt.cell <- 0
                
                }
                
                # Define probabilities for each model
                prob.wt <- c(0.999, 0.5, 0.001)
                prob.mt <- c(0.001, 0.5, 0.999)
                
                scores <- NULL
                
                # Test each model
                for(i in 1:length(prob.wt)) {
                
                    # Retrieve test statistic
                    x2 <- chisq.test(c(counts.wt.cell, counts.mt.cell), p=c(prob.wt[i], prob.mt[i]))$statistic
                    
                    # Convert to zygosity score and save into vector
                    scores[i] <- 1/log10(x2 + 1)
                    
                }
                
                # Save into list
                scores.list[[j]] <- scores
                
                #print(j)
                    
            }

            # Tabulate results
            results <- do.call(rbind.data.frame, scores.list)

            # Assign genotype
            . <- apply(results, 1, function(x) {which.max(x)})
            .[which(.== "1")] <- "WT"
            .[which(.== "2")] <- "Het"
            .[which(.== "3")] <- "Hom"
            df.small..$genotype <- .
            
            # Save into list
            .list[[k]] <- df.small..
        
        } else {
        
            # Save into list
            .list[[k]] <- NULL
        
        }
        
        
    }

    results <- do.call(rbind.data.frame, .list)
    print(nrow(results)==nrow(df.small.))
    .list.del[[1]] <- results
    
}
    
# Assign genotype: VCF records + low coverage
df.small. <- df.small[which(!is.na(df.small$chr)), ]
#df.small. <- df.small.[which(df.small.$coverage <= df.small.$coverage.blank), ]
df.small. <- df.small.[which(df.small.$coverage <= 10), ]

if(nrow(df.small.) != 0) {

    .list <- list()

    variant_names <- unique(df.small.$variant_name)
    
    for(k in 1:length(variant_names)) {

        # Subset variant
        df.small.. <- df.small.[which(df.small.$variant_name==variant_names[k]), ]
        
        if(nrow(df.small..) != 0) {
        
            # Retrieve ref, alt nt
            nt.ref <- unique(md.variant[which(md.variant$variant_name==variant_names[k]), "alt"])
            nt.alt <- "del"
                
            # Retrieve ref, alt counts
            count.ref <- df.small..[, nt.ref]
            count.alt <- df.small..[, nt.alt]
            
            # Compute VAF
            df.small..$vaf <- count.alt / (count.ref + count.alt) * 100
            
            # Save into list
            .list[[k]] <- df.small..

        } else {
        
            # Save into list
            .list[[k]] <- NULL
        
        }

    }
    
    results <- do.call(rbind.data.frame, .list)
    print(nrow(results)==nrow(df.small.))
    df.small. <- results
    df.small.$genotype <- "Low Coverage"
    .list.del[[2]] <- df.small.
    
}
    
# Assign genotype: No VCF records
df.small. <- df.small[which(is.na(df.small$chr)), ]

if(nrow(df.small.) != 0) {
    
    df.small.$vaf <- 0
    #df.small.$genotype <- ifelse(df.small.$coverage > df.small.$coverage.blank, "WT", "Low Coverage")
    df.small.$genotype <- ifelse(df.small.$coverage > 10, "WT", "Low Coverage")
    .list.del[[3]] <- df.small.

}

# Merge
#sum(sapply(.list.del, function(x) {nrow(x)})) == nrow(df.small)
df.del <- do.call(rbind.data.frame, .list.del)

############################################################
####################### MERGE ##############################
############################################################

# Merge
nrow(df.ins) + nrow(df.del) == nrow(df)

df <- rbind.data.frame(df.ins, df.del)

# Indicate ambiguous genotype
table(df$genotype)
df$genotype[which(df$vaf > 2 & df$vaf < 4 & df$genotype != "Low Coverage")] <- "Ambiguous"
df$genotype[which(df$vaf > 96 & df$vaf < 98 & df$genotype != "Low Coverage")] <- "Ambiguous"
table(df$genotype)

# Recode missing VAF
#df$vaf[is.na(df$vaf)] <- 1

# Save file
path <- "/Users/seanwen/Documents/Vladimir/out/"
file <- "GST010_Mutect2_AlleleCounts_GenotypeAssigned_gDNA.txt"
write.table(df, paste(path, file, sep=""), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
