# Load Necessary Libraries
library(TxDb.Hsapiens.UCSC.hg38.knownGene)  # UCSC known gene annotation for Homo sapiens
library(readr)
library(tximport)
library(data.table)
library(tidyverse)
library(biomaRt)      

# Initialize variables for human gene annotations
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
tx2gene <- na.omit(tx2gene)  # Remove any rows with missing values
tx2gene  # View the resulting tx2gene DataFrame

# Load Meta Data
coldata_csv <- read_csv("/PATH/TO/INPUT/coldata.csv")
#> colnames(coldata_csv)
#[1] "files"     "names"     "condition"

coldata_mTOR <- data.frame(files=coldata_csv$files, 
                           names=coldata_csv$names, 
                           condition=coldata_csv$condition, 
                           stringsAsFactors=FALSE)

# Create a vector of file paths
files_mTOR <- c(coldata_mTOR$files)
all(file.exists(files_mTOR)) # Ensure all files exist

# --------------------------------------------------------------------
# Optional for Salmon: Clean first column of Salmon output to extract ENST only; easier for tximport
# Define a function to modify each file
# If you do not want to overwrite files modify the function
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
txi_mTOR <- tximport(files_mTOR, type = "salmon", tx2gene = tx2gene)
txi_mTOR.tx <- tximport(files_mTOR, type = "salmon", txOut = TRUE)
txi_mTOR.sum <- summarizeToGene(txi_mTOR.tx, tx2gene)
all.equal(txi_mTOR$counts, txi_mTOR.sum$counts)

# Write the abundance data to a DataFrame
read_abundance_mTOR <- data.frame(txi_mTOR$abundance)
oldcolnames <- colnames(read_abundance_mTOR)
newcolnames <- coldata_mTOR$condition
setnames(read_abundance_mTOR, old = oldcolnames, new = newcolnames)
read_abundance_mTOR <- rownames_to_column(read_abundance_mTOR, var = "entrezgene_id")

# Map IDs using biomaRt
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
id_map <- getBM(
  attributes=c('ensembl_gene_id', 'external_gene_name', 'entrezgene_id'),
  filters = 'entrezgene_id',
  values = read_abundance_mTOR$entrezgene_id,
  mart = ensembl
)

# Merge mapped IDs with abundance data
read_abundance_mTOR <- read_abundance_mTOR %>%  mutate(entrezgene_id = as.numeric(entrezgene_id))
read_abundance_mTOR <- left_join(id_map, read_abundance_mTOR, by = "entrezgene_id")

# Rename external_gene_name to Symbol
read_abundance_mTOR <- read_abundance_mTOR %>% rename("Symbol" = "external_gene_name")

# Remove duplicate symbols and set columns as numeric
read_abundance_mTOR <- read_abundance_mTOR[!duplicated(read_abundance_mTOR$ensembl_gene_id), ]
read_abundance_mTOR <- read_abundance_mTOR[!duplicated(read_abundance_mTOR$Symbol), ]
read_abundance_mTOR[newcolnames] <- sapply(read_abundance_mTOR[newcolnames],as.numeric)

# Save the abundance data to a CSV file
write_csv(read_abundance_mTOR, "/PATH/TO/OUTPUT/mTOR_abundance.csv")

# Import filtered non-coding RNA with BASH (“MIR-” “SNOR-” “LINC-” and “LNC-”)
read_abundance_mTOR <- read.csv("/PATH/TO/INPUT/mTOR_abundance_filtered.csv")
TableOfCounts_mTOR <- read_abundance_mTOR[-c(1,3)]
write_csv(TableOfCounts_mTOR, "/PATH/TO/OUTPUT/TableOfCounts_mTOR_filtered.csv")

# Modify table of counts for further analysis
TableOfCounts_mTOR_modified <- data.frame(TableOfCounts_mTOR, row.names = 1)

# Check for any missing values
sum(is.na(TableOfCounts_mTOR_modified))
