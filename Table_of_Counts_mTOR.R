# Load Necessary Libraries
library(biomaRt)      
library(TxDb.Hsapiens.UCSC.hg38.knownGene)  # UCSC known gene annotation for Homo sapiens
library(readr)        
library(dplyr)
library(tximport)
library(tibble)
library(stringr)

# Initialize variables for human gene annotations
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
tx2gene <- na.omit(tx2gene)  # Remove any rows with missing values
tx2gene  # View the resulting tx2gene DataFrame

# Load Meta Data
coldata_csv <- read_csv("coldata.csv")

coldata_mtor <- data.frame(files=coldata_csv$files, 
                           names=coldata_csv$names, 
                           condition=coldata_csv$condition, 
                           stringsAsFactors=FALSE)

# Create a vector of file paths
files_mtor <- c(coldata_mtor$files)
all(file.exists(files_mtor)) # Ensure all files exist

# --------------------------------------------------------------------
# Optional: Clean first column of Salmon output to extract ENST only
# Define a function to modify each file
modify_file <- function(file_path) {
  df <- read.table(file_path, header = TRUE, sep = "\t")
  df$Name <- sub("^(ENST\\d+\\.\\d+).*", "\\1", df$Name)  # Extract ENST IDs
  return(df)
}

# Loop through each file, modify it, and save (overwrite) the modified version
for (file_path in files) {
  modified_df <- modify_file(file_path)
  write.table(modified_df, file_path, sep = "\t", row.names = FALSE, quote = FALSE)
}
# --------------------------------------------------------------------

# Import Transcript Quantification Data with tximport
txi_mtor <- tximport(files_mtor, type = "salmon", tx2gene = tx2gene)
txi_mtor.tx <- tximport(files_mtor, type = "salmon", txOut = TRUE)
txi_mtor.sum <- summarizeToGene(txi_mtor.tx, tx2gene)
all.equal(txi_mtor$counts, txi_mtor.sum$counts)

# Write the abundance data to a DataFrame
read_abundance_mtor <- data.frame(txi_mtor$abundance)
oldcolnames <- colnames(read_abundance_mtor)
newcolnames <- coldata_mtor$condition
setnames(read_abundance_mtor, old = oldcolnames, new = newcolnames)
read_abundance_mtor <- rownames_to_column(read_abundance_mtor, var = "entrezgene_id")

# Map IDs using biomaRt
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
id_map <- getBM(
  attributes=c('ensembl_gene_id', 'external_gene_name', 'entrezgene_id'),
  filters = 'entrezgene_id',
  values = read_abundance_atf4$entrezgene_id,
  mart = ensembl
)

# Merge mapped IDs with abundance data
read_abundance_mtor <- read_abundance_mtor %>%  mutate(entrezgene_id = as.numeric(entrezgene_id))
read_abundance_mtor <- left_join(id_map, read_abundance_mtor, by = "entrezgene_id")

# Rename external_gene_name to Symbol
read_abundance_mtor <- read_abundance_mtor %>% rename("external_gene_name" = "Symbol")

# Remove duplicate symbols and set columns as numeric
read_abundance_mtor <- read_abundance_mtor[!duplicated(read_abundance_mtor$ensembl_gene_id), ]
read_abundance_mtor <- read_abundance_mtor[!duplicated(read_abundance_mtor$Symbol), ]
read_abundance_mtor[newcolnames] <- sapply(read_abundance_mtor[newcolnames],as.numeric)

# Save the abundance data to a CSV file
write_csv(read_abundance_mtor, "mTOR_abundance.csv")

# Import filtered non-coding RNA except SNRP with BASH (“MIR-” “SNOR-” “LINC-” and “LNC-”)
read_abundance_mtor <- read.csv("/PATH/TO/INPUT/mTOR_abundance_filtered.csv")
TableOfCounts_mtor <- read_abundance_mtor[-c(1,3)]
write_csv(TableOfCounts_mtor, "TableOfCounts_mTOR_filtered.csv")

# Modify table of counts for further analysis
TableOfCounts_mtor_modified <- data.frame(TableOfCounts_mtor, row.names = 1)

# Check for any missing values
sum(is.na(TableOfCounts_mtor_modified))
