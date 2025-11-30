#!/usr/bin/env Rscript

# Create transcript-to-gene mapping from GTF file or cDNA FASTA headers
# Usage: Rscript create_tx2gene.R <gtf_or_fasta_file> <output_file>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript create_tx2gene.R <gtf_or_fasta_file> <output_file>")
}

input_file <- args[1]
output_file <- args[2]

cat("Creating tx2gene mapping...\n")
cat("Input file:", input_file, "\n")
cat("Output file:", output_file, "\n\n")

# Check if input is FASTA (easier and more reliable for Ensembl cDNA)
if (grepl("\\.fa$|\\.fasta$|\\.fa\\.gz$|\\.fasta\\.gz$", input_file, ignore.case = TRUE)) {
  cat("Detected FASTA file - parsing headers...\n")
  
  # Read FASTA file
  if (grepl("\\.gz$", input_file)) {
    lines <- readLines(gzfile(input_file))
  } else {
    lines <- readLines(input_file)
  }
  
  # Get header lines
  headers <- lines[grepl("^>", lines)]
  
  # Parse Ensembl cDNA headers
  # Format: >ENST00000456328.2 cdna chromosome:GRCh38:1:11869:14409:1 gene:ENSG00000223972.5 ...
  tx2gene_list <- lapply(headers, function(h) {
    # Extract transcript ID (after >)
    tx_match <- regmatches(h, regexpr("^>([^ ]+)", h, perl = TRUE))
    transcript_id <- sub("^>", "", tx_match)
    # Remove version number if present
    transcript_id <- sub("\\.[0-9]+$", "", transcript_id)
    
    # Extract gene ID
    gene_match <- regmatches(h, regexpr("gene:([^ ]+)", h, perl = TRUE))
    if (length(gene_match) > 0) {
      gene_id <- sub("gene:", "", gene_match)
      # Remove version number if present
      gene_id <- sub("\\.[0-9]+$", "", gene_id)
      return(c(transcript_id, gene_id))
    } else {
      return(NULL)
    }
  })
  
  # Remove NULL entries and convert to data frame
  tx2gene_list <- tx2gene_list[!sapply(tx2gene_list, is.null)]
  tx2gene <- as.data.frame(do.call(rbind, tx2gene_list), stringsAsFactors = FALSE)
  colnames(tx2gene) <- c("transcript_id", "gene_id")
  
} else {
  # Assume GTF file
  cat("Detected GTF file - parsing annotations...\n")
  
  # Read GTF file
  gtf <- read.table(input_file, sep = "\t", comment.char = "#", 
                    stringsAsFactors = FALSE, quote = "")
  colnames(gtf) <- c("seqname", "source", "feature", "start", "end", 
                     "score", "strand", "frame", "attributes")
  
  # Filter for transcript features
  gtf_tx <- gtf[gtf$feature == "transcript", ]
  
  if (nrow(gtf_tx) == 0) {
    cat("Warning: No transcript features found. Trying 'exon' features instead...\n")
    gtf_tx <- gtf[gtf$feature == "exon", ]
  }
  
  # Extract transcript_id and gene_id from attributes using multiple patterns
  extract_attribute <- function(attr_string, attr_name) {
    # Try with quotes
    pattern1 <- paste0(attr_name, ' "([^"]+)"')
    matches <- regmatches(attr_string, regexpr(pattern1, attr_string, perl = TRUE))
    if (length(matches) > 0 && matches != "") {
      return(sub(paste0('.*', attr_name, ' "([^"]+)".*'), "\\1", attr_string))
    }
    # Try without quotes
    pattern2 <- paste0(attr_name, ' ([^;]+)')
    matches <- regmatches(attr_string, regexpr(pattern2, attr_string, perl = TRUE))
    if (length(matches) > 0 && matches != "") {
      return(sub(paste0('.*', attr_name, ' ([^;]+).*'), "\\1", attr_string))
    }
    return(NA)
  }
  
  transcript_ids <- sapply(gtf_tx$attributes, extract_attribute, "transcript_id")
  gene_ids <- sapply(gtf_tx$attributes, extract_attribute, "gene_id")
  
  # Remove version numbers
  transcript_ids <- sub("\\.[0-9]+$", "", transcript_ids)
  gene_ids <- sub("\\.[0-9]+$", "", gene_ids)
  
  # Create tx2gene data frame
  tx2gene <- data.frame(
    transcript_id = transcript_ids,
    gene_id = gene_ids,
    stringsAsFactors = FALSE
  )
}

# Remove NA values
tx2gene <- tx2gene[complete.cases(tx2gene), ]

# Remove duplicates
tx2gene <- unique(tx2gene)

cat("Created", nrow(tx2gene), "transcript-to-gene mappings\n")

if (nrow(tx2gene) == 0) {
  stop("ERROR: No mappings created. Please check your input file format.")
}

# Preview first few mappings
cat("\nFirst few mappings:\n")
print(head(tx2gene, 5))

# Write output
write.table(tx2gene, output_file, sep = "\t", quote = FALSE, 
            row.names = FALSE, col.names = FALSE)

cat("\nMapping file created successfully!\n")
