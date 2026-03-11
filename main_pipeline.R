# ==============================================================================
# PROJECT: TCGA MULTI-OMICS INTEGRATION (Final Corrected)
# GOAL:    Identify Basal vs Luminal A Signatures using DIABLO
# AUTHOR:  Mohamed Sayed Ahmed
# ==============================================================================

# --- MODULE 0: ENVIRONMENT SETUP ---
rm(list = ls())      # Clean workspace memory
graphics.off()       # Close any stuck plot windows

# 1. Define Professional Output Directories
base_path  <- file.path(Sys.getenv("USERPROFILE"), "Documents")
data_dir   <- file.path(base_path, "TCGA_DATA_CACHE")      
output_dir <- file.path(base_path, "FINAL_ANALYSIS_FIGURES") 

if(!dir.exists(data_dir))   dir.create(data_dir, recursive = TRUE)
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 2. Load Core Bioinformatics Libraries
required_libs <- c("curatedTCGAData", "MultiAssayExperiment", "TCGAbiolinks", 
                   "mixOmics", "ggplot2", "reshape2", "ggpubr")

for(pkg in required_libs){
  if(!require(pkg, character.only = TRUE, quietly = TRUE)){
    if(pkg %in% c("ggplot2", "reshape2", "ggpubr")) install.packages(pkg)
    else BiocManager::install(pkg, update = FALSE)
  }
}

# --- MODULE 1: DATA INGESTION (SMART CACHING) ---
message(">>> [MODULE 1] Fetching TCGA-BRCA Data...")

local_file <- file.path(data_dir, "TCGA_BRCA_Raw.RData")

if(file.exists(local_file)) {
  message("Local cache found! Loading data from disk (Instant)...")
  load(local_file)
} else {
  message("No local cache found. Downloading from Cloud (This happens only once)...")
  mae <- curatedTCGAData(diseaseCode = "BRCA", assays = c("RNASeq2GeneNorm", "miRNASeqGene"), 
                         dry.run = FALSE, version = "2.0.1")
  save(mae, file = local_file)
  message("Data downloaded and saved to: ", local_file)
}

# --- MODULE 2: DATA ALIGNMENT & CLEANING ---
message(">>> [MODULE 2] Aligning Patients & Filtering Noise...")

# 1. Extract Matrices
exps <- experiments(mae)
ids_rna <- substr(colnames(exps[[1]]), 1, 12)
ids_mir <- substr(colnames(exps[[2]]), 1, 12)
common_patients <- intersect(ids_rna, ids_mir)

# 2. Fetch Clinical Subtypes
subtypes_db <- PanCancerAtlas_subtypes()
clinical_brca <- subtypes_db[subtypes_db$cancer.type == "BRCA", ]
clinical_brca$shortID <- substr(clinical_brca$pan.samplesID, 1, 12)

# 3. Filter for Target Analysis Groups
matched_clinical <- clinical_brca[clinical_brca$shortID %in% common_patients, ]
target_metadata  <- matched_clinical[matched_clinical$Subtype_mRNA %in% c("Basal", "LumA"), ]
final_ids        <- target_metadata$shortID

# 4. Remove Duplicates
if(any(duplicated(final_ids))) {
  message("Warning: Duplicates found. Removing matches...")
  keep_unique <- !duplicated(final_ids)
  final_ids <- final_ids[keep_unique]
  target_metadata <- target_metadata[keep_unique, ]
}

# 5. Create Aligned Matrices
X_mrna  <- assays(exps[[1]])[[1]][, match(final_ids, ids_rna)]
X_mirna <- assays(exps[[2]])[[1]][, match(final_ids, ids_mir)]

# Force unique row names
colnames(X_mrna) <- final_ids
colnames(X_mirna) <- final_ids

# --- [CRITICAL FIX] VARIANCE FILTERING ---
message("   ... Removing features with zero variance...")

# Function to identify non-informative genes (Variance = 0)
keep_var_rna <- apply(X_mrna, 1, var) > 0
keep_var_mir <- apply(X_mirna, 1, var) > 0

# Apply Filter
X_mrna <- X_mrna[keep_var_rna, ]
X_mirna <- X_mirna[keep_var_mir, ]

message("Cleaning Complete.")
message("   - Final Patients: ", length(final_ids))
message("   - Final Genes: ", nrow(X_mrna))
message("   - Final miRNAs: ", nrow(X_mirna))

# --- MODULE 3: MULTI-OMICS INTEGRATION (DIABLO) ---
message(">>> [MODULE 3] Running DIABLO Integration...")

# Prepare Input List (Transpose so Rows = Patients)
X <- list(mRNA = t(X_mrna), miRNA = t(X_mirna))
Y <- as.factor(target_metadata$Subtype_mRNA)

# Define Design Matrix
design <- matrix(0.1, ncol = length(X), nrow = length(X), dimnames = list(names(X), names(X)))
diag(design) <- 0

# Train DIABLO Model (Now safe from Zero Variance errors)
MyResult.diablo <- block.splsda(X, Y, ncomp = 2, design = design,
                                keepX = list(mRNA = c(20, 20), miRNA = c(20, 20)))

message("Model Trained Successfully!")

# --- MODULE 4: VISUALIZATION ---
message(">>> [MODULE 4] Generating Publication Figures...")
pub_colors <- c("firebrick3", "steelblue4")

# A. CIRCOS PLOT
tiff(file.path(output_dir, "Figure1_Circos.tiff"), width=10, height=10, units="in", res=300)
circosPlot(MyResult.diablo, cutoff = 0.7, color.blocks = pub_colors, size.labels = 1.5)
dev.off()

# B. SCALED INTEGRATED HEATMAP
top_genes  <- selectVar(MyResult.diablo, comp = 1)$mRNA$name
top_mirnas <- selectVar(MyResult.diablo, comp = 1)$miRNA$name
combined_mat <- scale(cbind(X$mRNA[, top_genes], X$miRNA[, top_mirnas]))

long_df <- melt(combined_mat)
colnames(long_df) <- c("Patient", "Feature", "ZScore")
long_df$Subtype <- Y[match(long_df$Patient, rownames(X$mRNA))]

p4 <- ggplot(long_df, aes(x = Patient, y = Feature, fill = ZScore)) +
  geom_tile() +
  facet_grid(~Subtype, scales = "free_x", space = "free") +
  scale_fill_gradient2(low = "navy", mid = "white", high = "firebrick3", limits = c(-2, 2), oob = scales::squish) +
  theme_minimal() + 
  theme(axis.text.x = element_blank(), panel.grid = element_blank()) +
  labs(title = "Multi-Omics Signature: Basal vs Luminal A", x = "Patients", y = "Integrated Features")

ggsave(file.path(output_dir, "Figure4_Integrated_Heatmap.pdf"), p4, width = 12, height = 8)

message("SUCCESS! All figures saved in: ", output_dir)
