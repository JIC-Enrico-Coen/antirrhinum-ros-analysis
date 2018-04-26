######################################
## Differential expression analysis ##
######################################
#
# Summarily, this script:
# - reads and compiles the transcript counts from SALMON software
# - summarises counts at the gene-level (using tximport)
# - exports table of normalised counts
# - exports a DESeq2 object 
# - Exports DESeq2 results for ROS/ros and EL/el contrasts

# load packages
library(DESeq2)
library(tidyverse)
library(tximport)

# set working directoy to project folder
##setwd("~/mount/hpc_slcu/jic/rosel_rnaseq/")
setwd("~/Documents/work/jic/hz_draft/supplementary_files/rnaseq_buds/")

#
# Data import ----
#
# List of all quantification files
quant_files <- list.files("data/counts/", pattern = "quant.sf", 
           full.names = TRUE,
           recursive = TRUE)
names(quant_files) <- dirname(quant_files) %>% basename()

# Create table with correspondence between gene and transcript names
## this is extracted from the name of the transcripts from one of the quant files
tx2gene <- quant_files[[1]] %>%
  read_tsv() %>%
  select(transcript = Name) %>%
  mutate(gene = str_replace(transcript, ".T.*$", ""))

# import with tximport
gene_quant <- tximport(quant_files, type = "salmon", tx2gene = tx2gene)

# Read sample information
sample_info <- read.csv("data/sample_info.csv", stringsAsFactors = FALSE) %>% 
  select(-fastq_folder) %>% 
  rename(library = tgac_library_name) %>% 
  arrange(library) %>%
  mutate(rosel_alleles = str_replace(rosel_alleles, " .*", ""))

# Ensure sample information and column names of counts are the same
if(all(colnames(gene_quant$counts) == sample_info$library)){
  # Name rows of sample_info
  rownames(sample_info) <- colnames(gene_quant$counts)
} else {
  stop("Samples do not match")
}


#
# Run DE tests ----
#
# Create DESeq2 object
rosel_dds <- DESeqDataSetFromTximport(gene_quant, colData = sample_info, design = ~ ros_genotype + el_genotype)

# Run test - we want to shrink log fold change estimates to de-emphasise genes with low counts (which are noisier)
rosel_dds <- DESeq(rosel_dds, betaPrior = TRUE)

# Produce variance-stabilised and normalised counts
rosel_rlog <- rlog(rosel_dds)


#
# Export full data ----
#
# Save the DESeq objects
write_rds(rosel_rlog, "./scripts/R/processed_data/DESeqTransform_all.rds", compress = "gz")
write_rds(rosel_dds, "./scripts/R/processed_data/DESeqDataSet_all.rds", "gz")

# TPMs
gene_quant$abundance %>% 
  as_tibble(rownames = "gene") %>% 
  gather(library, cts, -gene) %>% 
  full_join(sample_info, by = "library") %>% 
  write_csv("./scripts/R/processed_data/counts_tpm_all.csv")

# variance stabilised normalised counts
assay(rosel_rlog) %>% 
  as_tibble(rownames = "gene") %>% 
  gather(library, cts, -gene) %>% 
  full_join(sample_info, by = "library") %>% 
  write_csv("./scripts/R/processed_data/counts_rlog_all.csv")

# DESeq results - save both ROS vs ros and EL vs el contrasts
list(c("ros_genotype", "ROS", "ros"), c("el_genotype", "EL", "el")) %>%  # These are the contrasts
  lapply(function(i){                                # loop through it to get all results
    
    # Compute and tidy results from Walt test
    results(rosel_dds, contrast = c(i[1], i[2], i[3])) %>% 
      as.data.frame() %>% as_tibble(rownames = "gene_id") %>% 
      mutate(contrast = paste(i[2], "vs", i[3]))
  
    }) %>% 
  bind_rows() %>% 
  mutate(padj2 = p.adjust(pvalue, method = "BH")) %>%     # make new adjusted p-value for all tests (more stringent)
  write_csv("./scripts/R/processed_data/deseq_tests_all.csv")



#
# Re-analysis without "outlier" sample ----
#
# The PCA reveals an outlier sample (one of the JI7 genotypes)
# We therefore will re-run the test excluding this sample
plotPCA(rosel_rlog, "rosel_alleles") + geom_label(aes(label = name))

# import with tximport
outlier_index <- which(names(quant_files) == "LIB6740")
gene_quant <- tximport(quant_files[-outlier_index], type = "salmon", tx2gene = tx2gene)

# Read sample information
sample_info <- sample_info %>% 
  filter(library != "LIB6740") %>% 
  droplevels()

# Ensure sample information and column names of counts are the same
if(all(colnames(gene_quant$counts) == sample_info$library)){
  # Name rows of sample_info
  rownames(sample_info) <- colnames(gene_quant$counts)
} else {
  stop("Samples do not match")
}


#
# Run DE tests ----
#
# Create DESeq2 object
rosel_dds <- DESeqDataSetFromTximport(gene_quant, colData = sample_info, design = ~ ros_genotype + el_genotype)

# Run test
rosel_dds <- DESeq(rosel_dds, betaPrior = TRUE)

# Produce variance-stabilised normalised counts
rosel_rlog <- rlog(rosel_dds)


#
# Export data excluding "outlier" sample ----
#
# Save the DESeq objects
write_rds(rosel_rlog, "./scripts/R/processed_data/DESeqTransform_outlier.rds", compress = "gz")
write_rds(rosel_dds, "./scripts/R/processed_data/DESeqDataSet_outlier.rds", "gz")

# TPMs
gene_quant$abundance %>% 
  as_tibble(rownames = "gene") %>% 
  gather(library, cts, -gene) %>% 
  full_join(sample_info, by = "library") %>% 
  write_csv("./scripts/R/processed_data/counts_tpm_outlier.csv")

# variance stabilised normalised counts
assay(rosel_rlog) %>% 
  as_tibble(rownames = "gene") %>% 
  gather(library, cts, -gene) %>% 
  full_join(sample_info, by = "library") %>% 
  write_csv("./scripts/R/processed_data/counts_rlog_outlier.csv")

# DESeq results - report all results
list(c("ros_genotype", "ROS", "ros"), c("el_genotype", "EL", "el")) %>%  # These are the contrasts
  lapply(function(i){                                # loop through it to get all results
    
    # Compute and tidy results from Walt test
    results(rosel_dds, contrast = c(i[1], i[2], i[3])) %>% 
      as.data.frame() %>% as_tibble(rownames = "gene_id") %>% 
      mutate(contrast = paste(i[2], "vs", i[3]))
    
  }) %>% 
  bind_rows() %>% 
  mutate(padj2 = p.adjust(pvalue, method = "BH")) %>%     # make new adjusted p-value for all tests
  write_csv("./scripts/R/processed_data/deseq_tests_outlier.csv")


