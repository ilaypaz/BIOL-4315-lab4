BiocManager::install("clusterProfiler")
BiocManager::install("EnhancedVolcano")
BiocManager::install("biomaRt", force = TRUE)
#Note, if biomaRt isn't installed on your machine for some reason, you should install it as well. 

lapply(c(
  "docopt", "DT", "pheatmap","GenomicFeatures", "DESeq2",
  "edgeR", "systemPipeR", "systemPipeRdata", "BiocStyle", "GO.db", "dplyr",
  "tidyr", "stringr", "Rqc", "QuasR", "DT", "ape", "clusterProfiler","biomaRt","EnhancedVolcano"
), require, character.only = TRUE)
library(biomaRt)
library(DT)

# Step 1: Get fastq path and file list (optional, for context)
fq_path <- systemPipeRdata::pathList()$fastqdir
fq_files <- list.files(fq_path)
print(fq_files)

# Step 2: Read the original metadata file (unaltered)
file_path <- "/Users/ipaz00/Library/R/arm64/4.4/library/systemPipeRdata/extdata/param/targetsPE.txt"
meta_data <- read.table(file_path, header = TRUE, stringsAsFactors = FALSE)
head(meta_data)

# Step 3: Simulate replacements IN MEMORY (no file is overwritten!)
file_content <- readLines(file_path)

# Replace *_1.fastq.gz → .fna path (just for display)
file_content_modified <- gsub(
  pattern = "\\./data/SRR446[0-9]+_1\\.fastq\\.gz",
  replacement = "/Users/ipaz00/Downloads/BIOL4315_R/biol4315lab4/BIOL-4315-lab4/GCF_000001735.4_TAIR10.1_genomic.fna",
  x = file_content
)

# Replace *_2.fastq.gz → .gtf path (just for display)
file_content_modified <- gsub(
  pattern = "\\./data/SRR446[0-9]+_2\\.fastq\\.gz",
  replacement = "/Users/ipaz00/Downloads/BIOL4315_R/biol4315lab4/BIOL-4315-lab4/GCF_000001735.4_TAIR10.1_genomic.gtf",
  x = file_content_modified
)

# Convert modified content into a data frame (skip comments if needed)
meta_data_modified <- read.table(
  text = file_content_modified, 
  header = TRUE, 
  stringsAsFactors = FALSE
)

# Step 4: Display in interactive table (modified version for viewing only)
DT::datatable(meta_data_modified)

#QUESTIon 2 


library(Rqc)
fastq_files <- unique(c(meta_data$FileName1, meta_data$FileName2))
qcRes <- rqc(path=fq_path, pattern = ".fastq.gz")
rqcCycleQualityBoxPlot(qcRes)
rqcReadFrequencyPlot(qcRes)


rqcCycleQualityBoxPlot(qcRes[0:12])
rqcCycleQualityBoxPlot(qcRes[13:24])
rqcCycleQualityBoxPlot(qcRes[25:36])

rqcReadFrequencyPlot(qcRes[0:12])
rqcReadFrequencyPlot(qcRes[13:24])
rqcReadFrequencyPlot(qcRes[25:36])

#QUESTION 3
library(ShortRead)
targets_data <- read.table(
  "/Users/ipaz00/Library/R/arm64/4.4/library/systemPipeRdata/extdata/param/targetsPE.txt",
  header = TRUE,stringsAsFactors = FALSE)
fastq_dir <- "/Users/ipaz00/Library/R/arm64/4.4/library/systemPipeRdata/extdata/fastq"
targets_data$FileName1 <- file.path(fastq_dir, basename(targets_data$FileName1))
targets_data$FileName2 <- file.path(fastq_dir, basename(targets_data$FileName2))
output_dir <- "processed_reads"
if (!dir.exists(output_dir)) dir.create(output_dir)
adapter_seq <- "GCCCGGGTAA"

process_read_file <- function(input_path, sample_name) {
  output_file <- file.path(output_dir, paste0(sample_name, "_processed.fastq.gz"))
  fq <- readFastq(input_path)
  trimmed_seq <- trimLRPatterns(Lpattern = adapter_seq, subject = sread(fq), max.Lmismatch = 0)
  keep_idx <- width(trimmed_seq) > 0
  fq <- fq[keep_idx]
  trimmed_seq <- trimmed_seq[keep_idx]
  trimmed_seq <- narrow(trimmed_seq, end = width(trimmed_seq) - 3)
  trimmed_qual <- narrow(quality(fq), end = width(trimmed_seq))  # Note: use trimmed_seq widths here to sync lengths
  fq <- ShortReadQ(sread = trimmed_seq, quality = trimmed_qual, id = ShortRead::id(fq)[keep_idx])
if(length(fq) == 0) {
  writeFastq(fq, output_file, compress = TRUE)
  return(output_file)
}
afreq <- alphabetFrequency(sread(fq), baseOnly = TRUE)
if(!"N" %in% colnames(afreq)) {
  keep_noN <- rep(TRUE, length(fq))
} else {
  keep_noN <- afreq[, "N"] <= 1
}
fq <- fq[keep_noN]

  writeFastq(fq, output_file, compress = TRUE)
  return(output_file)
}

processed_files <- mapply(
  FUN = process_read_file,
  input_path = targets_data$FileName1,
  sample_name = targets_data$SampleName,
  SIMPLIFY = TRUE
)

targets_data$ProcessedFile1 <- processed_files

DT::datatable(targets_data)

#whatever this is 
#create the output directory
dir.create("./hisat2_index", recursive = TRUE)

at_genome <- "GCF_000001735.4_TAIR10.1_genomic.fna"
#Use system2 to run hisat2 from within R

tryCatch({system2(command = "hisat2-build", 
                  args = c("-p","8", at_genome,
                           "./hisat2_index/tair10_1_index"),
                  stdout = TRUE, stderr = TRUE)}, error = function(e) {
                    paste("hisat2-build", "indexing failed with error:", e$message)
                  })

sam_dir <- "/Users/ipaz00/Downloads/BIOL4315_R/biol4315lab4/BIOL-4315-lab4/sam_files"
bam_dir <- "/Users/ipaz00/Downloads/BIOL4315_R/biol4315lab4/BIOL-4315-lab4/bam_files"
fastq_dir <- "/Users/ipaz00/Downloads/BIOL4315_R/biol4315lab4/BIOL-4315-lab4/processed_reads"
hisat2_index <- "/Users/ipaz00/Downloads/BIOL4315_R/biol4315lab4/BIOL-4315-lab4/hisat2_index/tair10_1_index"


meta_data <- read.table("/Users/ipaz00/Library/R/arm64/4.4/library/systemPipeRdata/extdata/param/targetsPE.txt", header = TRUE, stringsAsFactors = FALSE)

for (sample_name in meta_data$SampleName) {
  fastq_file <- file.path(fastq_dir, paste0(sample_name, "_processed.fastq"))
  sam_file <- file.path(sam_dir, paste0(sample_name, ".sam"))
  
  if (!file.exists(fastq_file)) {
    message("Skipping ", sample_name, ": FASTQ file not found - ", fastq_file)
    next
  }
  
  hisat2_log <- tryCatch({
    system2("hisat2", args = c("-x", hisat2_index, "-U", fastq_file, "-S", sam_file), stdout = TRUE, stderr = TRUE)
  }, error = function(e) paste("hisat2 failed:", e$message))
  writeLines(hisat2_log, file.path(bam_dir, paste0(sample_name, "_hisat2.log")))
  
  bam_file <- file.path(bam_dir, paste0(sample_name, ".bam"))
  sorted_bam_file <- file.path(bam_dir, paste0(sample_name, "_sorted.bam"))
  
  system2("samtools", args = c("view", "-bS", sam_file, "-o", bam_file))
  system2("samtools", args = c("sort", "-o", sorted_bam_file, bam_file))
  system2("samtools", args = c("index", sorted_bam_file))
  
  cat("Finished processing sample:", sample_name, "\n")
}

library(dplyr)
library(stringr)

bam_dir <- "/Users/ipaz00/Downloads/BIOL4315_R/biol4315lab4/BIOL-4315-lab4/bam_files"
log_files <- list.files(bam_dir, pattern = "_hisat2\\.log$", full.names = TRUE)

if(length(log_files) == 0) stop("No HISAT2 log files found in bam_files directory.")

percent_aligned <- vector("character", length(log_files))

for (i in seq_along(log_files)) {
  lines <- readLines(log_files[i])
  matched_line <- grep("overall alignment rate", lines, value = TRUE)
  if(length(matched_line) == 1) {
    percent_aligned[i] <- matched_line
  } else {
    percent_aligned[i] <- NA
  }
}

align_df <- data.frame(
  samplename = sort(meta_data$SampleName),
  percent_aligned = percent_aligned,
  stringsAsFactors = FALSE
)

align_df <- align_df %>%
  mutate(
    percent_aligned = as.numeric(str_extract(percent_aligned, "\\d+\\.?\\d*"))
  )

head(align_df)

library(dplyr)
library(stringr)
library(DT)
library(ggplot2)
bam_dir <- "/Users/ipaz00/Downloads/BIOL4315_R/biol4315lab4/BIOL-4315-lab4/bam_files"
log_files <- list.files(bam_dir, pattern = "_hisat2\\.log$", full.names = TRUE)
log_sample_names <- basename(log_files) %>%
str_remove("_hisat2\\.log$") %>%
sort()
percent_aligned <- vector("character", length(log_files))
for (i in seq_along(log_files)) {
  lines <- readLines(log_files[i])
  matched_line <- grep("overall alignment rate", lines, value = TRUE)
  if(length(matched_line) == 1) {
    percent_aligned[i] <- matched_line
  } else {
    percent_aligned[i] <- NA
  }
}
align_df <- data.frame(
  samplename = log_sample_names,
  percent_aligned = percent_aligned,
  stringsAsFactors = FALSE)

align_df <- align_df %>%
  mutate(
    percent_aligned = as.numeric(str_extract(percent_aligned, "\\d+\\.?\\d*")))

head(align_df)
datatable(align_df)
ggplot(align_df, aes(x = "", y = percent_aligned)) +
  geom_boxplot(fill = "red") +
  labs(title =  "HISAT2 Alignment Percentages",
       y = "Alignment Percent (%)",x = "") 
##post question 5 pain 
BiocManager::install("Rsubread")
library(Rsubread)

# Define the full path to your bam files directory
bam_dir <- "/Users/ipaz00/Downloads/BIOL4315_R/biol4315lab4/BIOL-4315-lab4/bam_files"

# List all sorted BAM files with full path
bfiles <- list.files(bam_dir, pattern = "_sorted.bam$", full.names = TRUE)

# Define path to annotation GTF file (update path as necessary)
gtf_file <- "/path/to/GCF_000001735.4_TAIR10.1_genomic.gtf"

# Run featureCounts
gene_count_list <- featureCounts(
  files = bfiles,
  annot.ext = gtf_file,
  isGTFAnnotationFile = TRUE,
  allowMultiOverlap = FALSE,
  isPairedEnd = FALSE,
  nthreads = 8,
  minMQS = 10,
  GTF.featureType = "exon",
  GTF.attrType = "gene_id"
)

# View summary of counts
head(gene_count_list$counts)









