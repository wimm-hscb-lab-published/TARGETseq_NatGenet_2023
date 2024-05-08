# Load modules
module load R-base/4.0.1
module load R-cbrg/202109
R

# Define folders/cell_ids
path <- "/project/meadlab/wwen/Vladimir/GST010/"
cell.ids <- list.files(path)

# Tabulate coverage
.list <- list()

pb <- txtProgressBar(1, length(cell.ids), style=3)

for(i in 1:length(cell.ids)) {

    # Define cell_id
    cell.id <- cell.ids[i]

    # Define file name
    path.in <- paste(path, cell.id, "/mRNA/", sep="")
    file.in <- paste(cell.id, "_MuTect2.vcf", sep="")
    
    # Read file
    df <- try(
          read.table(paste(path.in, file.in, sep=""), sep="\t", header=FALSE, stringsAsFactors=FALSE),
          TRUE)
        
    if(class(df) != "try-error") {
            
        # Annotate cell id
        . <- data.frame("cell.id"=cell.id, stringsAsFactors=FALSE)
        df <- cbind.data.frame(., df)
        
        # Save into list
        .list[[i]] <- df
    
    }
                    
    # Track progress
    setTxtProgressBar(pb, i)
    
}

df <- do.call(rbind.data.frame, .list)

# Generate column names
names(df) <- c("cell.id", "chr", "start", "end", "ref", "alt", "strand", "score", "info", "format", "value")

# Correct mis-labelled T
df$ref[which(df$ref=="TRUE")] <- "T"
df$alt[which(df$alt=="TRUE")] <- "T"

# Save file
path <- "/project/meadlab/wwen/Vladimir/out/"
file <- "GST010_Mutect2_mRNA.txt"
write.table(df, paste(path, file, sep=""), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
