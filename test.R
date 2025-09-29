BiocManager::install("clusterProfiler")
BiocManager::install("EnhancedVolcano")
BiocManager::install("biomaRt", force = TRUE)
BiocManager::install("clusterProfiler", force=TRUE)
BiocManager::install("ShortRead", force = TRUE)
BiocManager::install("Rsubread",force=TRUE)
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
rqcCycleBaseCallsPlot(qcRes[0:12]) + theme_bw()
rqcCycleBaseCallsPlot(qcRes[13:24]) + theme_bw()
rqcCycleBaseCallsPlot(qcRes[25:36]) + theme_bw()
rqcReadFrequencyPlot(qcRes[0:12])
rqcReadFrequencyPlot(qcRes[13:24])
rqcReadFrequencyPlot(qcRes[25:36])


file.exists("processed_reads/M1A_processed.fastq.gz")
file.remove("processed_reads/M1A_processed.fastq.gz")
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
  trimmed_qual <- narrow(quality(fq), end = width(trimmed_seq))
  fq <- ShortReadQ(sread = trimmed_seq, quality = trimmed_qual, id = ShortRead::id(fq)[keep_idx])
  
  if(length(fq) == 0) {
    if(file.exists(output_file)) file.remove(output_file)
    writeFastq(fq, output_file, compress = TRUE, mode="w")
    return(output_file)
  }
  
  afreq <- alphabetFrequency(sread(fq), baseOnly = TRUE)
  if(!"N" %in% colnames(afreq)) {
    keep_noN <- rep(TRUE, length(fq))
  } else {
    keep_noN <- afreq[, "N"] <= 1
  }
  fq <- fq[keep_noN]
  
  if(file.exists(output_file)) file.remove(output_file)   # <-- add this line to avoid error
  writeFastq(fq, output_file, compress = TRUE, mode="w")
  return(output_file)
}


processed_reads <- mapply(
  FUN = process_read_file,
  input_path = targets_data$FileName1,
  sample_name = targets_data$SampleName,
  SIMPLIFY = TRUE
)

targets_data$ProcessedFile1 <- processed_reads

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
#question 4
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
library(Rsubread)

# Define the full path to your bam files directory
bam_dir <- "/Users/ipaz00/Downloads/BIOL4315_R/biol4315lab4/BIOL-4315-lab4/bam_files"

# List all sorted BAM files with full path
bfiles <- list.files(bam_dir, pattern = "_sorted.bam$", full.names = TRUE)

# Define path to annotation GTF file (update path as necessary)
gtf_file <- "/Users/ipaz00/Downloads/BIOL4315_R/biol4315lab4/BIOL-4315-lab4/GCF_000001735.4_TAIR10.1_genomic.gtf"

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

#list bam files
bfiles <- list.files("./bam_files", pattern = "_sorted.bam$", full.names = TRUE)

#Counting how many reads corespond to each gene
gene_count_list <- Rsubread::featureCounts(files = bfiles, annot.ext = "GCF_000001735.4_TAIR10.1_genomic.gtf", 
                                           isGTFAnnotationFile = TRUE, # <--- input annotation is GTF
                                           allowMultiOverlap = FALSE,  # <--- don't allow reads that overlap with multiple loci
                                           isPairedEnd = FALSE, nthreads = 8,
                                           minMQS = 10, # <--- minimum mapping quality score of 10 (like a phred score for the hisat2 alignment)
                                           GTF.featureType = "exon",  # <--- Count reads overlapping 'exon' features
                                           GTF.attrType = "gene_id" # <--- Groups exon by gene id
)
library(tibble)
tibble::glimpse(gene_count_list$annotation)[1:5,]


count_table <- gene_count_list$count

colnames(count_table) <- gsub("_sorted\\.bam$", "", colnames(count_table))


count_table <- count_table[rowSums(count_table[ , -1]) > 0, ]

datatable(count_table, 
              options = list(pageLength = 10), 
              caption = "Filtered Count Table")


#Get the metadata table that would accompany the count table
coldata <- meta_data %>% dplyr::select(SampleName,SampleLong,Factor) %>% 
  dplyr::mutate(SampleLong=str_split_i(SampleLong, "\\.",1)) %>% #getting the groups name (Avirulent, Mock and Virulent)
  dplyr::rename(condition = SampleLong) %>%
  dplyr::mutate(condition = factor(condition)) %>% #for group comp
  dplyr::mutate(Factor = factor(Factor)) #for sample comp

base::rownames(coldata) <- coldata$SampleName
coldata <- coldata %>% mutate(SampleName = factor(SampleName))

coldata$type <- factor(rep("single-read", nrow(coldata)))

#Make sure samples are in the same order between the two. 
coldata <- coldata[base::match(base::colnames(count_table), rownames(coldata)),]

#Check that they are indeed same order
all(rownames(coldata) == base::colnames(count_table))

#creating a dds object where the conditions are avr vs mock vs vir

dds1 <- DESeqDataSetFromMatrix(countData = count_table,
                               colData = coldata,
                               design = ~ condition)
print(dds1)
#creating a dds object where the conditions are the samples themslevs 
dds2 <- DESeqDataSetFromMatrix(countData = count_table,
                               colData = coldata,
                               design = ~ Factor)
print(dds2)

#correlating the samples
d <- cor(assay(rlog(dds1)), method = "spearman")
#turning correlation to a distance (1 - correlation) and clustring
hc <- hclust(dist(1 - d))

#hirarchal clustering
plot.phylo(as.phylo(hc), type = "p", edge.col = "blue", edge.width = 2,
           show.node.label = TRUE, no.margin = TRUE)
dds1_results <- DESeq(dds1)
dds2_results <- DESeq(dds2)
res1 <- DESeq2::results(dds1_results)
res1
res2 <- DESeq2::results(dds2_results)
res2
dds2
res_vir_mock <- DESeq2::results(dds1_results, contrast = c("condition", "Vir", "Mock"), alpha = 0.2)

# Avr vs Mock  
res_avr_mock <- DESeq2::results(dds1_results, contrast = c("condition", "Avr", "Mock"), alpha = 0.2)

# Vir vs Avr
res_vir_avr <- DESeq2::results(dds1_results, contrast = c("condition", "Vir", "Avr"), alpha = 0.2)

# Function to filter and count DE genes
filter_and_count <- function(res_obj, comparison_name, fc_threshold = 2) {
  # Remove NAs
  res_filtered <- res_obj[!is.na(res_obj$padj) & !is.na(res_obj$log2FoldChange), ]
  
  # Apply filters: |log2FC| >= log2(2) = 1 and padj <= alpha (already set in results())
  sig_genes <- res_filtered[abs(res_filtered$log2FoldChange) >= log2(fc_threshold), ]
  
  # Count up and down regulated
  up_regulated <- sum(sig_genes$log2FoldChange > 0)
  down_regulated <- sum(sig_genes$log2FoldChange < 0)
  
  return(data.frame(
    Comparison = comparison_name,
    Up_regulated = up_regulated,
    Down_regulated = down_regulated
  ))
}

# Apply filtering and counting to all comparisons
results_summary <- rbind(
  filter_and_count(res_vir_mock, "Vir vs Mock"),
  filter_and_count(res_avr_mock, "Avr vs Mock"),
  filter_and_count(res_vir_avr, "Vir vs Avr")
)

# Print summary
print("Summary of DE genes (FC >= 2, alpha = 0.2):")

print(results_summary)

plot_data <- results_summary %>%
  pivot_longer(cols = c(Up_regulated, Down_regulated), 
               names_to = "Regulation", 
               values_to = "Count") %>%
  mutate(Regulation = factor(Regulation, levels = c("Up_regulated", "Down_regulated")))

# Create horizontal stacked bar plot
p <- ggplot(plot_data, aes(x = Comparison, y = Count, fill = Regulation)) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +  # Makes it horizontal
  labs(
    title = "Differentially Expressed Genes by Comparison",
    subtitle = "Fold Change >= 2, alpha = 0.2",
    x = "Comparison",
    y = "Number of Genes",
    fill = "Regulation"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10))

#QUESTION 7 
comp <- systemPipeR::readComp("/Users/ipaz00/Library/R/arm64/4.4/library/systemPipeRdata/extdata/param/targetsPE.txt")
pairs<-comp[[1]]
pairs
results_list <- list()
for (pair in pairs) {
  groups <- unlist(strsplit(pair, "-"))
  g1 <- groups[1]
  g2 <- groups[2]
  comp_name <- paste(g1, "vs", g2)
  res <- DESeq2::results(dds2_results, contrast = c("Factor", g1, g2), alpha = 0.2)
  filtered_summary <- filter_and_count(res, comp_name)
  results_list[[comp_name]] <- filtered_summary
}

pairwise_summary <- do.call(rbind, results_list)

plot_data <- pairwise_summary %>%
  pivot_longer(cols = c(Up_regulated, Down_regulated),
               names_to = "Regulation",
               values_to = "Count") %>%
  mutate(Regulation = factor(Regulation, levels = c("Up_regulated", "Down_regulated")))

p2 <- ggplot(plot_data, aes(x = Comparison, y = Count, fill = Regulation)) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +
  labs(
    title = "Differentially Expressed Genes by Sample Comparison",
    subtitle = "Using dds2 (Factor-level), FC ≥ 2, alpha = 0.2",
    x = "Comparison",
    y = "Number of DE Genes",
    fill = "Regulation"
  ) +
  theme_minimal()
theme(
  plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
  plot.subtitle = element_text(hjust = 0.5, size = 12),
  axis.text = element_text(size = 10),
  axis.title = element_text(size = 12),
  legend.title = element_text(size = 11),
  legend.text = element_text(size = 10))

p2

m <- biomaRt::useMart("plants_mart", dataset = "athaliana_eg_gene",
                      host = "https://plants.ensembl.org")
desc <- biomaRt::getBM(attributes = c("tair_locus", "description"), mart = m)
desc <- desc[!duplicated(desc[, 1]), ]
desc <- desc %>% rename( gene_id = tair_locus)

library(stringr)

annotate_results <- function(res_obj, desc_df) {
  res_df <- as.data.frame(res_obj)
  res_df$gene_id <- rownames(res_df)
  res_df <- left_join(res_df, desc_df, by = "gene_id") %>%
    mutate(description = str_replace(description, "\\[.*", "")) # removes everything from [ onward
  
  return(res_df)
}

# Add annotations to all results
res_vir_mock_annot <- annotate_results(res_vir_mock, desc)
res_avr_mock_annot <- annotate_results(res_avr_mock, desc)
res_vir_avr_annot <- annotate_results(res_vir_avr, desc)

res_vir_mock_annot$description[is.na(res_vir_mock_annot$description)] <- ""
res_vir_mock_annot <- res_vir_mock_annot %>% 
  filter(!is.na(pvalue) & !is.na(log2FoldChange))
library(dplyr)

# Replace NA in description with empty string
res_vir_mock_annot$description[is.na(res_vir_mock_annot$description)] <- ""

# Filter out NA in pvalue and log2FoldChange
res_vir_mock_annot <- res_vir_mock_annot %>%
  filter(!is.na(pvalue) & !is.na(log2FoldChange))

# Create short_label with first 25 characters of description
res_vir_mock_annot <- res_vir_mock_annot %>%
  mutate(short_label = substr(description, 1, 25)) %>%
  mutate(label = ifelse(abs(log2FoldChange) > 1 & pvalue < 0.05, short_label, ""))

# Plot with boxedLabels = FALSE to start, to avoid viewport errors
volcano1 <- EnhancedVolcano(res_vir_mock_annot,
                            lab = res_vir_mock_annot$label,
                            x = 'log2FoldChange',
                            y = 'pvalue',
                            title = 'Vir vs Mock',
                            pCutoff = 0.05,
                            FCcutoff = 1.0,
                            pointSize = 4.0,
                            labSize = 4.0,
                            labCol = 'black',
                            labFace = 'bold',
                            boxedLabels = FALSE,       # start with FALSE to avoid layout issues
                            colAlpha = 4/5,
                            legendPosition = 'right',
                            legendLabSize = 14,
                            legendIconSize = 4.0,
                            drawConnectors = TRUE,
                            widthConnectors = 1.0,
                            colConnectors = 'black')

print(volcano1)


# Step 1: Replace NA descriptions with empty strings
res_vir_avr_annot$description[is.na(res_vir_avr_annot$description)] <- ""

# Step 2: Create truncated labels (first 25 characters)
res_vir_avr_annot$short_label <- substr(res_vir_avr_annot$description, 1, 25)

# Step 3: Plot with simplified options first
volcano2_fixed <- EnhancedVolcano(res_vir_avr_annot,
                                  lab = res_vir_avr_annot$short_label,
                                  x = 'log2FoldChange',
                                  y = 'pvalue',
                                  title = 'Vir vs avir',
                                  pCutoff = 0.05,
                                  FCcutoff = 1.0,
                                  pointSize = 4.0,
                                  labSize = 4.0,
                                  labCol = 'black',
                                  labFace = 'bold',
                                  boxedLabels = TRUE,        # disable boxed labels first
                                  colAlpha = 4/5,
                                  legendPosition = 'right',
                                  legendLabSize = 14,
                                  legendIconSize = 4.0,
                                  drawConnectors = TRUE,     # disable connectors first
                                  widthConnectors = 1.0,
                                  colConnectors = 'black') +
  ggplot2::scale_y_continuous(breaks = seq(0,7,1))

print(volcano2_fixed)

# Replace NAs in description
res_avr_mock_annot$description[is.na(res_avr_mock_annot$description)] <- ""

# Create short labels
res_avr_mock_annot$short_label <- substr(res_avr_mock_annot$description, 1, 25)

# Step 1: Basic plot with simplified labeling options
volcano3_fixed <- EnhancedVolcano(res_avr_mock_annot,
                                  lab = res_avr_mock_annot$short_label,
                                  x = 'log2FoldChange',
                                  y = 'pvalue',
                                  title = 'Avr vs Mock',
                                  pCutoff = 0.05,
                                  FCcutoff = 1.0,
                                  pointSize = 4.0,
                                  labSize = 4.0,
                                  labCol = 'black',
                                  labFace = 'bold',
                                  boxedLabels = TRUE,    # disable boxed labels
                                  colAlpha = 4/5,
                                  legendPosition = 'right',
                                  legendLabSize = 14,
                                  legendIconSize = 4.0,
                                  drawConnectors = TRUE,  # disable connectors
                                  widthConnectors = 1.0,
                                  colConnectors = 'black') +
  ggplot2::scale_y_continuous(breaks = seq(0,4,1))

print(volcano3_fixed)








