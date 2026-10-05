# ==========================
# Load required packages
# ==========================
library(fgsea)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(tibble)

# ==========================
# Set working directory
# ==========================
# Replace with the directory containing the input files
dir <- "path/to/neuron"

# Files containing the complete PLS-ranked gene lists
ASD_files <- c("ASD1", "ASD2")

# Functional gene-set files
category_files <- c(
  "Axon",
  "bume",
  "carbo",
  "cell-pro",
  "Dendrite",
  "differentiation",
  "migra",
  "mye",
  "pota",
  "synapse"
)

# Output directory
out_dir <- file.path(dir, "fgsea_results_15632")
dir.create(
  out_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


# ==========================
# 1. Load AHBA background genes
# ==========================

# The file should contain 15,632 AHBA genes in the first column
AHBA_genes <- readLines(
  file.path(dir, "AHBA_15632.txt"),
  encoding = "UTF-8"
)

# Remove the header
AHBA_genes <- AHBA_genes[AHBA_genes != "GENE"]

# Remove leading/trailing spaces and empty lines
AHBA_genes <- trimws(AHBA_genes)
AHBA_genes <- AHBA_genes[AHBA_genes != ""]

# Remove duplicated genes
AHBA_genes <- unique(AHBA_genes)

length(AHBA_genes)

cat(
  "Number of AHBA background genes:",
  length(AHBA_genes),
  "\n"
)

if (length(AHBA_genes) != 15632) {
  warning(
    paste0(
      "AHBA background contains ",
      length(AHBA_genes),
      " genes instead of 15632."
    )
  )
}


# ==========================
# 2. Load functional gene sets
# ==========================
category_list <- lapply(category_files, function(f) {
  
  genes <- read.table(
    file.path(dir, f),
    header = TRUE,
    stringsAsFactors = FALSE
  )[, 1]
  
  # Remove duplicated genes
  genes <- unique(genes)
  
  # Retain only genes included in the 15,632 AHBA genes
  genes <- genes[genes %in% AHBA_genes]
  
  return(genes)
})

names(category_list) <- category_files


# Check the number of genes from each gene set
# that are represented in the AHBA background
gene_set_sizes <- data.frame(
  pathway = names(category_list),
  setSize_in_AHBA = sapply(category_list, length)
)

write.csv(
  gene_set_sizes,
  file.path(
    out_dir,
    "gene_set_sizes_in_AHBA15632.csv"
  ),
  row.names = FALSE
)


# ==========================
# 3. Define FGSEA analysis function
# ==========================
run_fgsea <- function(file_name) {
  
  cat("\n============================\n")
  cat("Running:", file_name, "\n")
  cat("============================\n")
  
  
  # --------------------------
  # Load the complete PLS results
  # --------------------------
  data <- read.table(
    file.path(dir, file_name),
    header = TRUE,
    stringsAsFactors = FALSE
  )
  
  # Expected columns:
  # GENE  = gene symbol
  # SCORE = PLS1 bootstrap Z-score
  ranks <- data$SCORE
  names(ranks) <- data$GENE
  
  
  # --------------------------
  # Remove missing values
  # --------------------------
  ranks <- ranks[!is.na(ranks)]
  ranks <- ranks[!is.na(names(ranks))]
  ranks <- ranks[names(ranks) != ""]
  
  
  # --------------------------
  # Retain only the 15,632
  # AHBA background genes
  # --------------------------
  ranks <- ranks[
    names(ranks) %in% AHBA_genes
  ]
  
  
  # --------------------------
  # Check duplicated genes
  # --------------------------
  if (anyDuplicated(names(ranks)) > 0) {
    
    cat(
      "Duplicated genes detected. ",
      "Keeping the first occurrence.\n"
    )
    
    ranks <- ranks[
      !duplicated(names(ranks))
    ]
  }
  
  
  # --------------------------
  # Check ranked-list size
  # --------------------------
  cat(
    "Number of genes in ranked list:",
    length(ranks),
    "\n"
  )
  
  if (length(ranks) != 15632) {
    
    warning(
      paste0(
        file_name,
        " contains ",
        length(ranks),
        " AHBA genes rather than 15632."
      )
    )
  }
  
  
  # --------------------------
  # Check missing AHBA genes
  # --------------------------
  missing_genes <- setdiff(
    AHBA_genes,
    names(ranks)
  )
  
  cat(
    "Missing AHBA genes:",
    length(missing_genes),
    "\n"
  )
  
  
  # --------------------------
  # Rank genes in descending order
  # --------------------------
  ranks <- sort(
    ranks,
    decreasing = TRUE
  )
  
  # ==========================
  # Run FGSEA
  # ==========================
  set.seed(12345)
  
  fgseaRes <- fgseaMultilevel(
    pathways = category_list,
    stats = ranks,
    minSize = 10,
    maxSize = 5000
  )
  
  
  # ==========================
  # Add leading-edge size
  # ==========================
  fgseaRes_out <- fgseaRes %>%
    mutate(
      
      leadingEdgeSize =
        sapply(
          leadingEdge,
          length
        ),
      
      leadingEdge =
        sapply(
          leadingEdge,
          function(x)
            paste(
              x,
              collapse = ","
            )
        )
    )
  
  
  # ==========================
  # Map pathway display names
  # ==========================
  pathway_name_map <- c(
    "Dendrite" = "Dendrite development",
    "bume" = "BSS",
    "carbo" = "Carbonic anhydrase",
    "cell-pro" = "Cell proliferation marker",
    "differentiation" = "Neuron differentiation",
    "migra" = "Neuron migration",
    "mye" = "Myelination",
    "pota" = "PCC",
    "synapse" = "Synapse development",
    "Axon" = "Axon development"
  )
  
  # Replace abbreviated names with display names
  fgseaRes_out$pathway <-
    pathway_name_map[
      fgseaRes_out$pathway
    ]
  
  
  # ==========================
  # Define pathway output order
  # ==========================
  pathway_order <- c(
    "Dendrite development",
    "BSS",
    "Carbonic anhydrase",
    "Cell proliferation marker",
    "Neuron differentiation",
    "Neuron migration",
    "Myelination",
    "PCC",
    "Synapse development",
    "Axon development"
  )
  
  fgseaRes_out <- fgseaRes_out %>%
    mutate(
      pathway = factor(
        pathway,
        levels = pathway_order
      )
    ) %>%
    arrange(pathway)
  
  
  # ==========================
  # Save FGSEA results
  # ==========================
  write.csv(
    fgseaRes_out,
    file = file.path(
      out_dir,
      paste0(
        file_name,
        "_fgsea_AHBA15632.csv"
      )
    ),
    row.names = FALSE
  )
  
  
  # ==========================
  # Generate NES heatmap
  # ==========================
  nes_mat <- fgseaRes_out %>%
    select(
      pathway,
      NES
    ) %>%
    mutate(
      pathway =
        as.character(pathway)
    ) %>%
    column_to_rownames(
      "pathway"
    )
  
  
  pdf(
    file.path(
      out_dir,
      paste0(
        file_name,
        "_NES_heatmap.pdf"
      )
    ),
    width = 5,
    height = 7
  )
  
  pheatmap(
    nes_mat,
    
    # Preserve the predefined pathway order
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    
    color =
      colorRampPalette(
        c(
          "blue",
          "white",
          "red"
        )
      )(50),
    
    main =
      paste0(
        file_name,
        " NES"
      )
  )
  
  dev.off()
  
  return(
    fgseaRes_out
  )
}


# ==========================
# 4. Run FGSEA for ASD1 and ASD2
# ==========================
fgsea_results <- lapply(
  ASD_files,
  run_fgsea
)

names(fgsea_results) <- ASD_files