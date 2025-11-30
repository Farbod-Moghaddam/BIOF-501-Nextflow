#!/usr/bin/env Rscript

# Merge Kallisto quantification results into gene-level count and TPM matrices
# Usage: Rscript merge_kallisto.R <tx2gene_file>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript merge_kallisto.R <tx2gene_file>")
}

tx2gene_file <- args[1]

cat("Merging Kallisto results...\n")
cat("tx2gene file:", tx2gene_file, "\n\n")

# Get all sample directories
sample_dirs <- list.dirs(recursive = FALSE, full.names = TRUE)
sample_names <- basename(sample_dirs)

cat("Found", length(sample_dirs), "sample directories:\n")
print(sample_names)
cat("\n")

# Load transcript to gene mapping
tx2gene <- read.table(tx2gene_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(tx2gene) <- c("transcript_id", "gene_id")
cat("Loaded", nrow(tx2gene), "transcript-to-gene mappings\n\n")

# Initialize lists to store data
all_counts <- list()
all_tpm <- list()

# Read each sample's abundance file
for (i in seq_along(sample_dirs)) {
    sample <- sample_names[i]
    abundance_file <- file.path(sample_dirs[i], "abundance.tsv")
    
    if (file.exists(abundance_file)) {
        cat("Processing", sample, "...\n")
        
        # Read kallisto output
        data <- read.table(abundance_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        cat("  Found", nrow(data), "transcripts\n")
        
        # Remove version numbers from transcript IDs (e.g., ENST00000123456.7 -> ENST00000123456)
        data$target_id_clean <- sub("\\..*", "", data$target_id)
        tx2gene$transcript_id_clean <- sub("\\..*", "", tx2gene$transcript_id)
        
        # Merge with gene mapping using cleaned IDs
        data_with_genes <- merge(data, tx2gene, by.x = "target_id_clean", by.y = "transcript_id_clean", all.x = TRUE)
        cat("  After merge:", nrow(data_with_genes), "rows\n")
        
        # Remove rows with missing gene_id
        data_with_genes <- data_with_genes[!is.na(data_with_genes$gene_id), ]
        cat("  After removing NA gene_ids:", nrow(data_with_genes), "rows\n")
        
        if (nrow(data_with_genes) > 0) {
            # Sum counts and TPM by gene (for transcripts mapping to same gene)
            gene_counts <- aggregate(est_counts ~ gene_id, data = data_with_genes, FUN = sum)
            gene_tpm <- aggregate(tpm ~ gene_id, data = data_with_genes, FUN = sum)
            
            # Store in lists
            all_counts[[sample]] <- gene_counts
            all_tpm[[sample]] <- gene_tpm
            cat("  Aggregated to", nrow(gene_counts), "genes\n")
        } else {
            stop("No valid gene mappings found for ", sample)
        }
    } else {
        stop("abundance.tsv not found for ", sample)
    }
}

cat("\nMerging all samples...\n")

# Merge all samples for counts
count_matrix <- all_counts[[1]]
colnames(count_matrix)[2] <- sample_names[1]
for (i in 2:length(all_counts)) {
    temp <- all_counts[[i]]
    colnames(temp)[2] <- sample_names[i]
    count_matrix <- merge(count_matrix, temp, by = "gene_id", all = TRUE)
}
count_matrix[is.na(count_matrix)] <- 0
rownames(count_matrix) <- count_matrix$gene_id
count_matrix <- count_matrix[, -1]

# Merge all samples for TPM
tpm_matrix <- all_tpm[[1]]
colnames(tpm_matrix)[2] <- sample_names[1]
for (i in 2:length(all_tpm)) {
    temp <- all_tpm[[i]]
    colnames(temp)[2] <- sample_names[i]
    tpm_matrix <- merge(tpm_matrix, temp, by = "gene_id", all = TRUE)
}
tpm_matrix[is.na(tpm_matrix)] <- 0
rownames(tpm_matrix) <- tpm_matrix$gene_id
tpm_matrix <- tpm_matrix[, -1]

# Save count matrix
write.table(count_matrix, "kallisto_count_matrix.txt", 
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
cat("Saved count matrix:", nrow(count_matrix), "genes x", ncol(count_matrix), "samples\n")

# Save TPM matrix
write.table(tpm_matrix, "kallisto_tpm_matrix.txt", 
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
cat("Saved TPM matrix:", nrow(tpm_matrix), "genes x", ncol(tpm_matrix), "samples\n")

cat("\nMerge completed successfully!\n")
