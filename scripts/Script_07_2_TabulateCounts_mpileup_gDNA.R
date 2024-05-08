# Load packages
library(data.table)
library(gdata)
library(plyr)
library(reshape2)

# Read files
    # VCF
    path <- "/Users/seanwen/Documents/Vladimir/out/"
    file <- "GST010_mpileup_gDNA.txt"
    df <- as.data.frame(fread(paste(path, file, sep=""), sep="\t", header=TRUE, stringsAsFactors=FALSE))
    
    # Variant metadata
    path <- "/Users/seanwen/Documents/Vladimir/Metadata/"
    file <- "Variants_Metadata.xlsx"
    md.variant <- read.xls(paste(path, file, sep=""), sheet=1, header=TRUE, stringsAsFactors=FALSE)
        
# Subset SNVs only
variant_types <- "SNV"
md.variant <- md.variant[which(md.variant$variant_type %in% variant_types), ]

############################ RECODE BASES #######################

bases.final <- NULL

pb <- txtProgressBar(1, nrow(df), style=3)

for(j in 1:nrow(df)) {

    # Retrieve bases
    bases <- df[j, "base"]
    bases <- toupper(bases)
        
    # Remove uninformative codings
        # Check full repotoire of codings
        # . <- paste(df$base, collapse="")
        # table(strsplit(., split=""))
        
        # Remove first base indicator
        bases <- gsub("\\^.", "", bases)
        
        # Remove last base indicator
        bases <- gsub("$", "", bases, fixed=TRUE)
        
        # Remove deleted base (after this read bases) indicators
        . <- strsplit(bases, split="")[[1]]
        
        del.logical <- grepl("-", ., fixed=TRUE)
        del.indices <- which(del.logical==TRUE)
        
        if(length(del.indices) != 0) {
        
            for(i in 1:length(del.indices)) {
                    
                del.logical <- grepl("-", ., fixed=TRUE)
                
                del.index <- which(del.logical==TRUE)[1]
                
                del.length.pos.1 <- as.numeric(.[del.index + 1])
                
                if(.[del.index + 2] == "N") {
                
                    del.length <- del.length.pos.1
                    
                    . <- .[-c(del.index, (del.index + 1):(del.index + 1 + del.length))]
                    
                } else {
                
                    del.length <- as.numeric(paste(.[c(del.index + 1):(del.index + 2)], collapse=""))
                    
                    . <- .[-c(del.index, (del.index + 1) ,(del.index + 2):(del.index + 2 + del.length))]
                
                
                }
                    
            
            
            }
            
            bases <- paste(., collapse="")
        
        }
        
        # Remove inserted base (after this read bases) indicators
        . <- strsplit(bases, split="")[[1]]
        
        ins.logical <- grepl("+", ., fixed=TRUE)
        ins.indices <- which(ins.logical==TRUE)
        
        if(length(ins.indices) != 0) {
        
            for(i in 1:length(ins.indices)) {
            
                ins.logical <- grepl("+", ., fixed=TRUE)
                
                ins.index <- which(ins.logical==TRUE)[1]
                
                ins.length <- as.numeric(.[ins.index + 1])
                
                . <- .[-c(ins.index, (ins.index + 1):(ins.index + 1 + ins.length))]
            
            }
            
            bases <- paste(., collapse="")
        
        }
        
    # Check final n bases matches coverage
    n.bases <- df[j, 3] == nchar(bases)
    
    if(n.bases==TRUE) {
    
        bases.final[j] <- bases
    
    
    } else {
    
        bases.final[j] <- bases
        print(j)
        
    }
    
    # Track progress
    setTxtProgressBar(pb, j)

}

# Trouble-shooting
# table(strsplit(bases, split=""))
# head(grep("-", df$base, fixed=TRUE))
# . <- paste(df$base, collapse="")
# table(strsplit(., split=""))

df$base <- bases.final

############################ COUTNS BY BASE TYPE #######################

# Check all base types
. <- paste(df$base, collapse="")
table(strsplit(., split=""))

# Tabulate counts
freq.list <- list()

pb <- txtProgressBar(1, nrow(df), style=3)

for(i in 1:nrow(df)) {

    # Tabulate frequency
    bases <- df[i, "base"]
    bases <- strsplit(bases, split="")
    freq <- as.data.frame(table(bases), stringsAsFactors=FALSE)
    
    # Recode deletion
    freq$bases[which(freq$bases=="*")] <- "del"
    
    # Create dummy entries for missing bases
    missing.bases <- setdiff(c("A", "C", "G", "T", "del", "ins"), unique(freq$bases))
    
    if(length(missing.bases)) {
    
    missing.freq <- data.frame("bases"=missing.bases,
                               "Freq"=0,
                               stringsAsFactors=FALSE
                               )
                               
    freq <- rbind.data.frame(freq, missing.freq)
    
    }
    
    # Set factor levels
    freq$bases <- factor(freq$bases, levels=c("A", "C", "G", "T", "del", "ins"))
    freq <- freq[order(freq$bases), ]
    
    # Reformat to single row
    freq <- as.data.frame(t(freq))
    names(freq) <- as.character(unlist(freq[1,]))
    freq <- freq[-1,]
    
    # Save into list
    freq.list[[i]] <- freq
    
    # Track progress
    setTxtProgressBar(pb, i)

}

freq.df <- do.call(rbind.data.frame, freq.list)
df <- cbind.data.frame(df, freq.df)
#df$base <- NULL
row.names(df) <- NULL

####################################################################################

# Annotate variant name
df <- join(df, md.variant[,c("start", "variant_name", "variant_type")], by="start", type="left")
sum(is.na(df$variant_name))

# Reorder columns to match Mutect2 output
df <- df[,c("cell.id", "chr", "start", "base", "A", "C", "G", "T", "ins", "del", "variant_name", "variant_type")]

# Save file
path <- "/Users/seanwen/Documents/Vladimir/out/"
file <- "GST010_mpileup_AlleleCounts_gDNA.txt"
write.table(df, paste(path, file, sep=""), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
