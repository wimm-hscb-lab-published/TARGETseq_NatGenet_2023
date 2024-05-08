# Load packages
library(data.table)
library(gdata)
library(plyr)
library(reshape2)

# Read files
    # VCF
    path <- "/Users/seanwen/Documents/Vladimir/out/"
    file <- "GST010_Mutect2_mRNA.txt"
    df <- as.data.frame(fread(paste(path, file, sep=""), sep="\t", header=TRUE, stringsAsFactors=FALSE))
    
    # Variant metadata
    path <- "/Users/seanwen/Documents/Vladimir/Metadata/"
    file <- "Variants_Metadata.xlsx"
    md.variant <- read.xls(paste(path, file, sep=""), sheet=1, header=TRUE, stringsAsFactors=FALSE)
        
# Subset INDELs only
variant_types <- c("Deletion", "Insertion")
md.variant <- md.variant[which(md.variant$variant_type %in% variant_types), ]

##########################

list.master <- list()

for(j in 1:nrow(md.variant)) {

    # Subset variant
    md.variant.small <- md.variant[j, ]
    
    index <- which(df$chr==md.variant.small$chr & df$start==md.variant.small$start)
    df.small <- df[index, ]
    
    if(nrow(df.small) != 0) {
    
        # Annotate variant type
            # Check nt change
            table(df.small$ref, df.small$alt)
            table(df.small$ref)
            table(df.small$alt)
            df.small$variant.type <- NA
            
            # Compute nt length
            df.small$ref.length <- as.character(nchar(df.small$ref))
            
            . <- strsplit(df.small$alt, split=",", fixed=TRUE)
            . <- sapply(., function(x) {nchar(x)})
            . <- sapply(., function(x) {paste(x, collapse="|")})
            df.small$alt.length <- as.character(.)
        
            table(df.small$ref.length, df.small$alt.length)
            
        # Assign: SNV
        . <- df.small[which(is.na(df.small$variant.type)), ]
        
        if(nrow(.) != 0) {
        
            logical <- NULL
            
            for(i in 1:nrow(.)) {
            
                ref <- .$ref[i]
                ref <- nchar(ref)
                
                alt <- .$alt[i]
                alt <- unlist(strsplit(alt, split=",", fixed=TRUE))
                alt <- nchar(alt)
                
                logical[[i]] <- length(which(alt==ref)) == length(alt)
                
            }
            
            cell.ids <- .$cell.id[logical==TRUE]
            
            df.small$variant.type[which(df.small$cell.id %in% cell.ids)] <- "SNV"
            
        }
            
        # Insertion
        . <- df.small[which(df.small$ref.length=="1" & is.na(df.small$variant.type)), ]
        
        if(nrow(.) != 0) {
        
            logical <- NULL
            
            for(i in 1:nrow(.)) {
            
                ref <- .$ref[i]
                ref <- nchar(ref)
                
                alt <- .$alt[i]
                alt <- unlist(strsplit(alt, split=",", fixed=TRUE))
                alt <- nchar(alt)
                
                if(length(unique(ref < alt)) != 1) {
                
                    logical[i] <- FALSE
                
                } else {
                
                    logical[i] <- unique(ref < alt) & length(unique(ref < alt))==1 & unique(ref < alt)==TRUE
                
                }
                
            }
            
            cell.ids <- .$cell.id[logical==TRUE]
            
            df.small$variant.type[which(df.small$cell.id %in% cell.ids)] <- "ins"
        
        }
            
        # Deletion
        . <- df.small[which(df.small$ref.length != "1" & is.na(df.small$variant.type)), ]
        
        if(nrow(.) != 0) {
        
            logical <- NULL
            
            for(i in 1:nrow(.)) {
            
                ref <- .$ref[i]
                ref <- nchar(ref)
                
                alt <- .$alt[i]
                alt <- unlist(strsplit(alt, split=",", fixed=TRUE))
                alt <- nchar(alt)
                
                if(length(unique(ref > alt)) != 1) {
                
                    logical[i] <- FALSE
                
                } else {
                
                    logical[i] <- unique(ref > alt) & length(unique(ref > alt))==1 & unique(ref > alt)==TRUE
                
                }
                
            }
        
        cell.ids <- .$cell.id[logical==TRUE]
        
        df.small$variant.type[which(df.small$cell.id %in% cell.ids)] <- "del"
        
        }
            
        # Assign: Complex: ins + del
        . <- df.small[which(is.na(df.small$variant.type)), ]
        
        if(nrow(.) != 0) {
        
            logical <- NULL
            
            for(i in 1:nrow(.)) {
            
                ref <- .$ref[i]
                ref <- nchar(ref)
                
                alt <- .$alt[i]
                alt <- unlist(strsplit(alt, split=",", fixed=TRUE))
                alt <- nchar(alt)
                
                logical[[i]] <- (length(which(alt < ref)) + length(which(alt > ref))) == length(alt)
                
            }
            
            cell.ids <- .$cell.id[logical==TRUE]
            
            df.small$variant.type[which(df.small$cell.id %in% cell.ids)] <- "Complex: ins + del"
            
        }
            
        # Assign: Complex: SNV + ins
        . <- df.small[which(is.na(df.small$variant.type)), ]
        
        if(nrow(.) != 0) {
        
            logical <- NULL
            
            for(i in 1:nrow(.)) {
            
                ref <- .$ref[i]
                ref <- nchar(ref)
                
                alt <- .$alt[i]
                alt <- unlist(strsplit(alt, split=",", fixed=TRUE))
                alt <- nchar(alt)
                
                logical[[i]] <- (length(which(alt==ref)) + length(which(alt > ref))) == length(alt)
                
            }
            
            cell.ids <- .$cell.id[logical==TRUE]
            
            df.small$variant.type[which(df.small$cell.id %in% cell.ids)] <- "Complex: SNV + ins"
            
        }
            
        # Assign: Complex: SNV + del
        . <- df.small[which(is.na(df.small$variant.type)), ]
        
        if(nrow(.) != 0) {
        
            logical <- NULL
            
            for(i in 1:nrow(.)) {
            
                ref <- .$ref[i]
                ref <- nchar(ref)
                
                alt <- .$alt[i]
                alt <- unlist(strsplit(alt, split=",", fixed=TRUE))
                alt <- nchar(alt)
                
                logical[[i]] <- (length(which(alt==ref)) + length(which(alt < ref))) == length(alt)
                
            }
            
            cell.ids <- .$cell.id[logical==TRUE]
            
            df.small$variant.type[which(df.small$cell.id %in% cell.ids)] <- "Complex: SNV + del"
            
        }
            
        # Assign: Complex: SNV + del
        . <- df.small[which(is.na(df.small$variant.type)), ]
        
        if(nrow(.) != 0) {
        
            logical <- NULL
            
            for(i in 1:nrow(.)) {
            
                ref <- .$ref[i]
                ref <- nchar(ref)
                
                alt <- .$alt[i]
                alt <- unlist(strsplit(alt, split=",", fixed=TRUE))
                alt <- nchar(alt)
                
                logical[[i]] <- (length(which(alt==ref)) + length(which(alt < ref)) + length(which(alt > ref))) == length(alt)
                
            }
            
            cell.ids <- .$cell.id[logical==TRUE]
            
            df.small$variant.type[which(df.small$cell.id %in% cell.ids)] <- "Complex: SNV + ins + del"
            
        }
        
        # Check for unannotated
        . <- df.small[is.na(df.small$variant.type), ]
        table(.$ref.length, .$alt.length)
        head(.)
        table(df.small$variant.type)
        sum(is.na(df.small$variant.type))
            
        # Sanity check
            # SNV
            . <- df.small[which(df.small$variant.type=="SNV"),]
            table(.$ref.length, .$alt.length)
            
            # Ins
            . <- df.small[which(df.small$variant.type=="ins"),]
            table(.$ref.length, .$alt.length)
            
            # Del
            . <- df.small[which(df.small$variant.type=="del"),]
            table(.$ref.length, .$alt.length)
            
            # Complex: ins + del
            . <- df.small[which(df.small$variant.type=="Complex: ins + del"),]
            table(.$ref.length, .$alt.length)
            
            # Complex: SNV + ins
            . <- df.small[which(df.small$variant.type=="Complex: SNV + ins"),]
            table(.$ref.length, .$alt.length)
            
            # Complex: SNV + del
            . <- df.small[which(df.small$variant.type=="Complex: SNV + del"),]
            table(.$ref.length, .$alt.length)
            
            # Complex: SNV + ins + del
            . <- df.small[which(df.small$variant.type=="Complex: SNV + ins + del"),]
            table(.$ref.length, .$alt.length)

        #################################### SNV ##########################################
        
        counts.list <- list()
        
        # Retrieve allele count
            # Subset relevant variant type
            df.small. <- df.small[which(df.small$variant.type=="SNV"), ]
            
            if(nrow(df.small.) != 0) {
            
                # Retrieve counts
                . <- strsplit(df.small.$value, split=":", fixed=TRUE)
                df.small.$count <- sapply(., function(x) {x[2]})
                        
                # Match counts w/ nt
                .list <- list()
                
                #pb <- txtProgressBar(1, nrow(df.small.), style=3)
                
                for(i in 1:nrow(df.small.)) {
                
                    # Retrieve individual counts
                    counts <- as.numeric(unlist(strsplit(df.small.$count[i], split=",", fixed=TRUE)))
                    
                    # Define ref allele
                    ref <- df.small.$ref[i]
                    nt.ref <- substr(ref, start=1, stop=1)
                    
                    # Define alt allele
                    alt <- df.small.$alt[i]
                    
                    if(length(grep(",", alt)) != 0) {
                    
                        alt <- unlist(strsplit(alt, split=",", fixed=TRUE))
                        nt.alt <- substr(alt, start=1, stop=1)
                        
                    }
                    
                    nt.alt <- substr(alt, start=1, stop=1)
                    
                    # Tabulate
                    .list[[i]] <- data.frame("cell.id"=df.small.$cell.id[i],
                                          "nt"=c(nt.ref, nt.alt),
                                          "count"=counts,
                                          stringsAsFactors=FALSE)
                                          
                    # Track progress
                    #setTxtProgressBar(pb, i)
                
                }
                
                results <- do.call(rbind.data.frame, .list)
                
                # Add dummy entry
                missing <- setdiff(c("A", "C", "G", "T", "ins", "del"), unique(results$nt))
                
                if(length(missing) != 0) {
                
                    missing.df <- data.frame("cell.id"="dummy",
                                             "nt"=missing,
                                             "count"=0,
                                             stringsAsFactors=FALSE)
                                             
                    results <- rbind.data.frame(results, missing.df)
                
                }
                
                # Set factor levels
                results$nt <- factor(results$nt, levels=c("A", "C", "G", "T", "ins", "del"))
                
                # Reshape data frame
                results <- dcast(data=results, formula=cell.id ~ nt, value.var="count", fun.aggregate=sum)
                
                # Annotate
                df.small. <- join(df.small., results, by="cell.id", type="left")
                
                # Save as new object
                df.small.snv <- df.small.
                
                # Save into list
                counts.list[[1]] <- df.small.
            
            }
        
        #################################### INS ##########################################
        
        # Retrieve allele count
            # Subset relevant variant type
            df.small. <- df.small[which(df.small$variant.type=="ins"), ]
            
            if(nrow(df.small.) != 0) {
            
                # Retrieve counts
                . <- strsplit(df.small.$value, split=":", fixed=TRUE)
                df.small.$count <- sapply(., function(x) {x[2]})
                        
                # Match counts w/ nt
                .list <- list()
                
                #pb <- txtProgressBar(1, nrow(df.small.), style=3)
                
                for(i in 1:nrow(df.small.)) {
                
                    # Retrieve individual counts
                    counts <- as.numeric(unlist(strsplit(df.small.$count[i], split=",", fixed=TRUE)))
                    
                    # Define ref allele
                    ref <- df.small.$ref[i]
                    nt.ref <- substr(ref, start=1, stop=1)
                    
                    # Define alt allele
                    alt <- df.small.$alt[i]
                    
                    n <- length(unlist(strsplit(alt, split=",", fixed=TRUE)))
                    
                    nt.alt <- rep("ins", times=n)
                            
                    # Tabulate
                    .list[[i]] <- data.frame("cell.id"=df.small.$cell.id[i],
                                          "nt"=c(nt.ref, nt.alt),
                                          "count"=counts,
                                          stringsAsFactors=FALSE)
                                          
                    # Track progress
                    #setTxtProgressBar(pb, i)
                    #print(i)
                    
                }
                
                results <- do.call(rbind.data.frame, .list)
                
                # Add dummy entry
                missing <- setdiff(c("A", "C", "G", "T", "ins", "del"), unique(results$nt))
                
                if(length(missing) != 0) {
                
                    missing.df <- data.frame("cell.id"="dummy",
                                             "nt"=missing,
                                             "count"=0,
                                             stringsAsFactors=FALSE)
                                             
                    results <- rbind.data.frame(results, missing.df)
                
                }
                
                # Set factor levels
                results$nt <- factor(results$nt, levels=c("A", "C", "G", "T", "ins", "del"))
                
                # Reshape data frame
                results <- dcast(data=results, formula=cell.id ~ nt, value.var="count", fun.aggregate=sum)
                
                # Annotate
                df.small. <- join(df.small., results, by="cell.id", type="left")
                
                # Save as new object
                df.small.ins <- df.small.
                
                # Save into list
                counts.list[[2]] <- df.small.
                
            }
        
        #################################### DEL ##########################################
        
        # Retrieve allele count
            # Subset relevant variant type
            df.small. <- df.small[which(df.small$variant.type=="del"), ]
            
            if(nrow(df.small.) != 0) {
            
                # Retrieve counts
                . <- strsplit(df.small.$value, split=":", fixed=TRUE)
                df.small.$count <- sapply(., function(x) {x[2]})
                        
                # Match counts w/ nt
                .list <- list()
                
                #pb <- txtProgressBar(1, nrow(df.small.), style=3)
                
                for(i in 1:nrow(df.small.)) {
                
                    # Retrieve individual counts
                    counts <- as.numeric(unlist(strsplit(df.small.$count[i], split=",", fixed=TRUE)))
                    
                    # Define ref allele
                    ref <- df.small.$ref[i]
                    nt.ref <- substr(ref, start=1, stop=1)
                    
                    # Define alt allele
                    alt <- df.small.$alt[i]
                    
                    n <- length(unlist(strsplit(alt, split=",", fixed=TRUE)))
                    
                    nt.alt <- rep("del", times=n)
                            
                    # Tabulate
                    .list[[i]] <- data.frame("cell.id"=df.small.$cell.id[i],
                                          "nt"=c(nt.ref, nt.alt),
                                          "count"=counts,
                                          stringsAsFactors=FALSE)
                                          
                    # Track progress
                    #setTxtProgressBar(pb, i)
                    #print(i)
                    
                }
                
                results <- do.call(rbind.data.frame, .list)
                
                # Add dummy entry
                missing <- setdiff(c("A", "C", "G", "T", "ins", "del"), unique(results$nt))
                
                if(length(missing) != 0) {
                
                    missing.df <- data.frame("cell.id"="dummy",
                                             "nt"=missing,
                                             "count"=0,
                                             stringsAsFactors=FALSE)
                                             
                    results <- rbind.data.frame(results, missing.df)
                
                }
                
                # Set factor levels
                results$nt <- factor(results$nt, levels=c("A", "C", "G", "T", "ins", "del"))
                
                # Reshape data frame
                results <- dcast(data=results, formula=cell.id ~ nt, value.var="count", fun.aggregate=sum)
                
                # Annotate
                df.small. <- join(df.small., results, by="cell.id", type="left")
                
                # Save as new object
                df.small.del <- df.small.
                
                # Save into list
                counts.list[[3]] <- df.small.
                
            }
        
        ########################### COMPLEX: INS + DEL #####################################
        
        # Retrieve allele count
            # Subset relevant variant type
            df.small. <- df.small[which(df.small$variant.type=="Complex: ins + del"), ]
            
            if(nrow(df.small.) != 0) {
            
                # Retrieve counts
                . <- strsplit(df.small.$value, split=":", fixed=TRUE)
                df.small.$count <- sapply(., function(x) {x[2]})
                        
                # Match counts w/ nt
                .list <- list()
                
                #pb <- txtProgressBar(1, nrow(df.small.), style=3)
                
                for(i in 1:nrow(df.small.)) {
                
                    # Retrieve individual counts
                    counts <- as.numeric(unlist(strsplit(df.small.$count[i], split=",", fixed=TRUE)))
                    
                    # Define ref allele
                    ref <- df.small.$ref[i]
                    nt.ref <- substr(ref, start=1, stop=1)
                    
                    # Define alt allele
                    alt <- df.small.$alt[i]
                    
                    alt <- unlist(strsplit(alt, split=",", fixed=TRUE))
                    
                    alt.ins.index <- which(nchar(ref) < nchar(alt))
                    alt.del.index <- which(nchar(ref) > nchar(alt))
                    
                    # Assign counts
                    counts.ref <- counts[1]
                    counts.alt.ins <- counts[alt.ins.index + 1]
                    counts.alt.del <- counts[alt.del.index + 1]
                    
                    nt.alt.ins <- rep("ins", times=length(counts.alt.ins))
                    nt.alt.del <- rep("del", times=length(counts.alt.del))
                            
                    # Tabulate
                    .list[[i]] <- data.frame("cell.id"=df.small.$cell.id[i],
                                          "nt"=c(nt.ref, nt.alt.ins, nt.alt.del),
                                          "count"=c(counts.ref, counts.alt.ins, counts.alt.del),
                                          stringsAsFactors=FALSE)
                                          
                    # Track progress
                    #setTxtProgressBar(pb, i)
                    #print(i)
                    
                }
                
                results <- do.call(rbind.data.frame, .list)
                
                # Add dummy entry
                missing <- setdiff(c("A", "C", "G", "T", "ins", "del"), unique(results$nt))
                
                if(length(missing) != 0) {
                
                    missing.df <- data.frame("cell.id"="dummy",
                                             "nt"=missing,
                                             "count"=0,
                                             stringsAsFactors=FALSE)
                                             
                    results <- rbind.data.frame(results, missing.df)
                
                }
                
                # Set factor levels
                results$nt <- factor(results$nt, levels=c("A", "C", "G", "T", "ins", "del"))
                
                # Reshape data frame
                results <- dcast(data=results, formula=cell.id ~ nt, value.var="count", fun.aggregate=sum)
                
                # Annotate
                df.small. <- join(df.small., results, by="cell.id", type="left")
                
                # Save as new object
                df.small.complex.ins.del <- df.small.
                
                # Save into list
                counts.list[[4]] <- df.small.
            
            }
        
        ########################### COMPLEX: SNV + INS #####################################
        
        # Retrieve allele count
            # Subset relevant variant type
            df.small. <- df.small[which(df.small$variant.type=="Complex: SNV + ins"), ]
            
            if(nrow(df.small.) != 0) {
            
                # Retrieve counts
                . <- strsplit(df.small.$value, split=":", fixed=TRUE)
                df.small.$count <- sapply(., function(x) {x[2]})
                        
                # Match counts w/ nt
                .list <- list()
                
                #pb <- txtProgressBar(1, nrow(df.small.), style=3)
                
                for(i in 1:nrow(df.small.)) {
                
                    # Retrieve individual counts
                    counts <- as.numeric(unlist(strsplit(df.small.$count[i], split=",", fixed=TRUE)))
                    
                    # Define ref allele
                    ref <- df.small.$ref[i]
                    nt.ref <- substr(ref, start=1, stop=1)
                    
                    # Define alt allele
                    alt <- df.small.$alt[i]
                    
                    alt <- unlist(strsplit(alt, split=",", fixed=TRUE))
                    
                    alt.snv.index <- which(nchar(ref) == nchar(alt))
                    alt.ins.index <- which(nchar(ref) < nchar(alt))
                    
                    # Assign counts
                    counts.ref <- counts[1]
                    counts.alt.snv <- counts[alt.snv.index + 1]
                    counts.alt.ins <- counts[alt.ins.index + 1]
                    
                    # Assign counts nt
                    nt.alt.snv <- substr(alt[alt.snv.index], start=1, stop=1)
                    nt.alt.ins <- rep("ins", times=length(counts.alt.ins))
                            
                    # Tabulate
                    .list[[i]] <- data.frame("cell.id"=df.small.$cell.id[i],
                                          "nt"=c(nt.ref, nt.alt.snv, nt.alt.ins),
                                          "count"=c(counts.ref, counts.alt.snv, counts.alt.ins),
                                          stringsAsFactors=FALSE)
                                          
                    # Track progress
                    #setTxtProgressBar(pb, i)
                    #print(i)
                    
                }
                
                results <- do.call(rbind.data.frame, .list)
                
                # Add dummy entry
                missing <- setdiff(c("A", "C", "G", "T", "ins", "del"), unique(results$nt))
                
                if(length(missing) != 0) {
                
                    missing.df <- data.frame("cell.id"="dummy",
                                             "nt"=missing,
                                             "count"=0,
                                             stringsAsFactors=FALSE)
                                             
                    results <- rbind.data.frame(results, missing.df)
                
                }
                
                # Set factor levels
                results$nt <- factor(results$nt, levels=c("A", "C", "G", "T", "ins", "del"))
                
                # Reshape data frame
                results <- dcast(data=results, formula=cell.id ~ nt, value.var="count", fun.aggregate=sum)
                
                # Annotate
                df.small. <- join(df.small., results, by="cell.id", type="left")
                
                # Save as new object
                df.small.complex.snv.ins <- df.small.

                # Save into list
                counts.list[[5]] <- df.small.
            
            }
        
        ########################### COMPLEX: SNV + DEL #####################################
        
        # Retrieve allele count
            # Subset relevant variant type
            df.small. <- df.small[which(df.small$variant.type=="Complex: SNV + del"), ]
            
            if(nrow(df.small.) != 0) {
            
                # Retrieve counts
                . <- strsplit(df.small.$value, split=":", fixed=TRUE)
                df.small.$count <- sapply(., function(x) {x[2]})
                        
                # Match counts w/ nt
                .list <- list()
                
                #pb <- txtProgressBar(1, nrow(df.small.), style=3)
                
                for(i in 1:nrow(df.small.)) {
                
                    # Retrieve individual counts
                    counts <- as.numeric(unlist(strsplit(df.small.$count[i], split=",", fixed=TRUE)))
                    
                    # Define ref allele
                    ref <- df.small.$ref[i]
                    nt.ref <- substr(ref, start=1, stop=1)
                    
                    # Define alt allele
                    alt <- df.small.$alt[i]
                    
                    alt <- unlist(strsplit(alt, split=",", fixed=TRUE))
                    
                    alt.snv.index <- which(nchar(ref) == nchar(alt))
                    alt.del.index <- which(nchar(ref) > nchar(alt))
                    
                    # Assign counts
                    counts.ref <- counts[1]
                    counts.alt.snv <- counts[alt.snv.index + 1]
                    counts.alt.del <- counts[alt.del.index + 1]
                    
                    # Assign counts nt
                    nt.alt.snv <- substr(alt[alt.snv.index], start=1, stop=1)
                    nt.alt.del <- rep("del", times=length(counts.alt.del))
                            
                    # Tabulate
                    .list[[i]] <- data.frame("cell.id"=df.small.$cell.id[i],
                                          "nt"=c(nt.ref, nt.alt.snv, nt.alt.del),
                                          "count"=c(counts.ref, counts.alt.snv, counts.alt.del),
                                          stringsAsFactors=FALSE)
                                          
                    # Track progress
                    #setTxtProgressBar(pb, i)
                    #print(i)
                    
                }
                
                results <- do.call(rbind.data.frame, .list)
                
                # Add dummy entry
                missing <- setdiff(c("A", "C", "G", "T", "ins", "del"), unique(results$nt))
                
                if(length(missing) != 0) {
                
                    missing.df <- data.frame("cell.id"="dummy",
                                             "nt"=missing,
                                             "count"=0,
                                             stringsAsFactors=FALSE)
                                             
                    results <- rbind.data.frame(results, missing.df)
                
                }
                
                # Set factor levels
                results$nt <- factor(results$nt, levels=c("A", "C", "G", "T", "ins", "del"))
                
                # Reshape data frame
                results <- dcast(data=results, formula=cell.id ~ nt, value.var="count", fun.aggregate=sum)
                
                # Annotate
                df.small. <- join(df.small., results, by="cell.id", type="left")
                
                # Save as new object
                df.small.complex.snv.del <- df.small.
                
                # Save into list
                counts.list[[6]] <- df.small.
            
            }
            
        ########################### COMPLEX: SNV + INS + DEL #####################################
        
        # Retrieve allele count
            # Subset relevant variant type
            df.small. <- df.small[which(df.small$variant.type=="Complex: SNV + ins + del"), ]
            
            if(nrow(df.small.) != 0) {
            
                # Retrieve counts
                . <- strsplit(df.small.$value, split=":", fixed=TRUE)
                df.small.$count <- sapply(., function(x) {x[2]})
                        
                # Match counts w/ nt
                .list <- list()
                
                #pb <- txtProgressBar(1, nrow(df.small.), style=3)
                
                for(i in 1:nrow(df.small.)) {
                
                    # Retrieve individual counts
                    counts <- as.numeric(unlist(strsplit(df.small.$count[i], split=",", fixed=TRUE)))
                    
                    # Define ref allele
                    ref <- df.small.$ref[i]
                    nt.ref <- substr(ref, start=1, stop=1)
                    
                    # Define alt allele
                    alt <- df.small.$alt[i]
                    
                    alt <- unlist(strsplit(alt, split=",", fixed=TRUE))
                    
                    alt.snv.index <- which(nchar(ref) == nchar(alt))
                    alt.ins.index <- which(nchar(ref) < nchar(alt))
                    alt.del.index <- which(nchar(ref) > nchar(alt))
                    
                    # Assign counts
                    counts.ref <- counts[1]
                    counts.alt.snv <- counts[alt.snv.index + 1]
                    counts.alt.ins <- counts[alt.ins.index + 1]
                    counts.alt.del <- counts[alt.del.index + 1]
                    
                    # Assign counts nt
                    nt.alt.snv <- substr(alt[alt.snv.index], start=1, stop=1)
                    nt.alt.ins <- rep("ins", times=length(counts.alt.ins))
                    nt.alt.del <- rep("del", times=length(counts.alt.del))
                            
                    # Tabulate
                    .list[[i]] <- data.frame("cell.id"=df.small.$cell.id[i],
                                          "nt"=c(nt.ref, nt.alt.snv, nt.alt.ins, nt.alt.del),
                                          "count"=c(counts.ref, counts.alt.snv, counts.alt.ins, counts.alt.del),
                                          stringsAsFactors=FALSE)
                                          
                    # Track progress
                    #setTxtProgressBar(pb, i)
                    #print(i)
                    
                }
                
                results <- do.call(rbind.data.frame, .list)
                
                # Add dummy entry
                missing <- setdiff(c("A", "C", "G", "T", "ins", "del"), unique(results$nt))
                
                if(length(missing) != 0) {
                
                    missing.df <- data.frame("cell.id"="dummy",
                                             "nt"=missing,
                                             "count"=0,
                                             stringsAsFactors=FALSE)
                                             
                    results <- rbind.data.frame(results, missing.df)
                
                }
                
                # Set factor levels
                results$nt <- factor(results$nt, levels=c("A", "C", "G", "T", "ins", "del"))
                
                # Reshape data frame
                results <- dcast(data=results, formula=cell.id ~ nt, value.var="count", fun.aggregate=sum)
                
                # Annotate
                df.small. <- join(df.small., results, by="cell.id", type="left")
                
                # Save as new object
                df.small.complex.snv.ins.del <- df.small.
                
                # Save into list
                counts.list[[7]] <- df.small.
            
            }
            
        #################################### MERGE ##########################################
        
        # Save final data frame
        results.final <- do.call(rbind.data.frame, counts.list)
        
        results.final$variant_name <- md.variant.small$variant_name
        results.final$variant_type <- md.variant.small$variant_type
                
        print(nrow(results.final)==nrow(df.small))
        
        list.master[[j]] <- results.final
        
        # Track progress
        print(paste(j, " of ", nrow(md.variant), " done", sep=""))
        
        
    } else {
    
        # Save final data frame
        list.master[[j]] <- NULL
        
        # Track progress
        print(paste(j, " of ", nrow(md.variant), " done", sep=""))
        
    }
        
        
}

# Merge all variants
sapply(list.master, function(x) {nrow(x)})
results <- do.call(rbind.data.frame, list.master)

# Save file
path <- "/Users/seanwen/Documents/Vladimir/out/"
file <- "GST010_Mutect2_AlleleCounts_mRNA.txt"
write.table(results, paste(path, file, sep=""), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
