# Load Necessary Libraries
library(biomaRt)      # BioMart data mining
library(TxDb.Hsapiens.UCSC.hg38.knownGene)  # UCSC known gene annotation for Homo sapiens
library(readr)        # Reading and writing CSV files
library(dplyr)        # Data manipulation
library(tximport)     # Import transcript quantification data
library(tibble)       # Data frame manipulation
library(stringr)      # String manipulation

# Initialize variables for human gene annotations
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
tx2gene <- na.omit(tx2gene)  # Remove any rows with missing values
tx2gene  # View the resulting tx2gene DataFrame

# Load Meta Data
coldata_csv <- read_csv("/Users/kwshin/Documents/capstone/ATF4KO/salmon_quant/coldata.csv")

# Convert metadata to a DataFrame
coldata_atf4 <- data.frame(
  files=coldata_csv$files,
  names=coldata_csv$names,
  condition=coldata_csv$condition,
  stringsAsFactors=FALSE
)

# Create a vector of file paths
files_atf4 <- c(coldata_atf4$files)
all(file.exists(files_atf4))  # Ensure all files exist

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
txi_atf4 <- tximport(files_atf4, type = "salmon", tx2gene = tx2gene)
txi_atf4_tx <- tximport(files_atf4, type = "salmon", txOut = TRUE)
txi_atf4_sum <- summarizeToGene(txi_atf4_tx, tx2gene)
all.equal(txi_atf4$counts, txi_atf4_sum$counts)  # Check equality of counts
attributes(txi_atf4)  # View attributes of the txi_atf4 object

# Write the abundance data to a DataFrame
read_abundance_atf4 <- data.frame(txi_atf4$abundance)
oldcolnames <- colnames(read_abundance_atf4)
newcolnames <- coldata_atf4$condition
setnames(read_abundance_atf4, old = oldcolnames, new = newcolnames)
read_abundance_atf4 <- rownames_to_column(read_abundance_atf4, var = "entrezgene_id")

# Map IDs using biomaRt
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
id_map <- getBM(
  attributes=c('ensembl_gene_id', 'external_gene_name', 'entrezgene_id'),
  filters = 'entrezgene_id',
  values = read_abundance_atf4$entrezgene_id,
  mart = ensembl
)

# Merge mapped IDs with abundance data
read_abundance_atf4 <- read_abundance_atf4 %>% mutate(entrezgene_id = as.numeric(entrezgene_id))
read_abundance_atf4 <- left_join(id_map, read_abundance_atf4, by = "entrezgene_id")

# Rename external_gene_name to Symbol
read_abundance_atf4 <- read_abundance_atf4 %>% rename("external_gene_name" = "Symbol")

# Remove duplicate symbols and set columns as numeric
read_abundance_atf4 <- read_abundance_atf4[!duplicated(read_abundance_atf4$ensembl_gene_id), ]
read_abundance_atf4 <- read_abundance_atf4[!duplicated(read_abundance_atf4$Symbol), ]
read_abundance_atf4[newcolnames] <- sapply(read_abundance_atf4[newcolnames], as.numeric)

# Save the cleaned abundance data to a CSV file
write_csv(read_abundance_atf4, "/PATH/TO/OUTPUT/ATF4KD_abundance.csv")

# Import filtered non-coding RNA except SNRP with BASH (“MIR-” “SNOR-” “LINC-” and “LNC-”)
read_abundance_atf4 <- read.csv("/PATH/TO/INPUT/ATF4KD_abundance_filtered.csv")
TableOfCounts_atf4 <- read_abundance_atf4[-c(1,3)]  # Remove first and third columns
write_csv(TableOfCounts_atf4, "/PATH/TO/OUTPUT/TableOfCounts_ATF4KD_filtered.csv")

# Modify table of counts for further analysis
TableOfCounts_atf4_modified <- data.frame(TableOfCounts_atf4, row.names = 1)

# Check for any missing values
sum(is.na(TableOfCounts_atf4_modified))
