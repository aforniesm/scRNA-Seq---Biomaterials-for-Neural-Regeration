############################################################
# 0. R PACKAGES NEEDED
############################################################

if (!require("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager")

# If the packages are not installed use this code:
# BiocManager::install(c("glmGamPoi", "biomaRt", "EnhancedVolcano", "speckle", "lisi", "AnnotationDbi","org.Mm.eg.db"))
# install.packages(c("Seurat", "Matrix", "tidyverse", "ggvenn", "tibble", "circlize"))

library(Seurat)
library(SeuratObject)
library(Matrix)
library(tidyverse)
library(glmGamPoi)    
library(EnhancedVolcano)
library(speckle)
library(lisi)
library(ggvenn)
library(ggplot2)
library(MAST)
library(ggrepel)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(biomaRt)
library(tibble)
library(circlize)

############################################################
# 1. LOAD THE DATA AND CREATE SEURAT OBJECT
############################################################

load_data <- function(matrix_path, 
                      metadata_path, 
                      features_path = NULL, 
                      barcodes_path = NULL, 
                      tech = "scRNA", 
                      study_name = "Study") {
  
  if (grepl("\\.mtx$", matrix_path)) {
    counts <- readMM(matrix_path)
    if (!is.null(features_path)) {
      features <- read.delim(features_path, header = FALSE)
      barcodes <- read.delim(barcodes_path, header = FALSE)
      rownames(counts) <- features[[1]]
      colnames(counts) <- barcodes[[1]]
    }
  } else if (grepl("\\.rds$", matrix_path)) {
    counts <- readRDS(matrix_path)
    if (class(counts) == "Seurat") counts <- GetAssayData(counts, "counts")
  } else {
    stop("Incorrect format for this Pipeline")
  }
  
  # Adjusting depending on the column names of the metadata
  metadata <- read.delim(metadata_path, header = TRUE) %>%
    rename(
      nCount_RNA = matches("nCount_RNA|nUMI"),
      nFeature_RNA = matches("nFeature_RNA|nGene"),
      celltype = matches("celltype|cluster|*cluster*"),
      time = matches("time|condition|label")
    ) %>%
  # Adjust depending on the conditions of the study and the different condition names
    mutate(
      time = case_when(
        grepl("(.*[^0-9]|^)7dpi|Acute|Contusion", time, ignore.case = TRUE) ~ "7dpi",
        grepl("(.*[^0-9]|^)1dpi", time, ignore.case = TRUE) ~ "1dpi",
        grepl("(.*[^0-9]|^)Uninjured|Control", time, ignore.case = TRUE) ~ "Uninjured",
        grepl("(.*[^0-9]|^)3.*wpi", time, ignore.case = TRUE) ~ "21dpi",
        grepl("(.*[^0-9]|^)6wpi|6wNT|SCI", time, ignore.case = TRUE) ~ "4-6wpi",
        TRUE ~ time
      ),
      technology = tech,
      dataset = study_name
    )
  
  # Create Seurat Object
  seurat <- CreateSeuratObject(
    counts = counts,
    meta.data = metadata,
    min.cells = 50,      
    min.features = 100  
  )
  
  # Quality control metrics
  removed_cells <- ncol(counts) - ncol(seurat)
  removed_features <- nrow(counts) - nrow(seurat)
  
  # Seurat object characteristics
  cells <- ncol(seurat)
  features <- ncol(seurat)
}

############################################################
# 2. NORMALIZATION
############################################################

normalize_data <- function(seurat_obj) {
  seurat_obj %>%
    NormalizeData(normalization.method = "LogNormalize") %>% 
    FindVariableFeatures(selection.method = "vst", nfeatures = 5000) 
}

############################################################
# 3. CELL TYPE ANNOTATION
############################################################

cell_annotation <- function(seurat_obj, reference) {
  anchors <- FindTransferAnchors(reference = reference, query = seurat_obj, dims = 1:20)
  predictions <- TransferData(anchorset = anchors, refdata = reference$celltype, dims = 1:20)
  seurat_obj <- AddMetaData(seurat_obj, metadata = predictions)
}

# Create a list with the Seurat Objects that you want to integrate
# seurat_list <- list()
  
############################################################
# 4. INTEGRATION
############################################################  

integration <- function(seurat_list, nfeatures) {
  features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = nfeatures)
  seurat_list <- lapply(X = seurat_list, FUN = function(x) {x <- ScaleData(x, features =  features)
  x <- RunPCA(x, features = features, npcs = 30)
  x <- RunUMAP(x, dims = 1:20)})
  anchors_integration <- FindIntegrationAnchors(object.list = seurat_list, dims = 1:20, verbose = TRUE, reduction = "rpca", anchor.features = features, scale = F, normalization.method = "LogNormalize")
  gc()
  seurat_integrated <- IntegrateData(anchorset = anchors_integration, dims = 1:20, verbose = TRUE, normalization.method = "LogNormalize")
  return(seurat_integrated)
}

############################################################
# *. UMAP BEFORE INTEGRATION (Add more maps if you integer more datasets)
############################################################  

umap_visualization <- function(seurat_list, path) {
  umap <- DimPlot(seurat_list [[1]], split.by = "dataset", group.by = "predicted.id", label = T, label.size = 3) + theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7)) + labs(x = "UMAP 1", y = "UMAP 2")
  umap2 <- DimPlot(seurat_list [[2]], split.by = "dataset", group.by = "predicted.id", label = T, label.size = 3) + theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7)) + labs(x = "UMAP 1", y = "UMAP 2")
  umap3 <- DimPlot(seurat_list [[3]], split.by = "dataset", group.by = "predicted.id", label = T, label.size = 3) + theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7)) + labs(x = "UMAP 1", y = "UMAP 2")
  umap + umap2 + umap3
  ggsave(path, umap + umap2 + umap3, dpi = 1080, height = 14, width = 8)
}

############################################################
# *. METADATA
############################################################  

metadata_cleaning <- function(seurat_obj) {
  metadata <- seurat_obj@meta.data
  metadata <- metadata %>% select(matches("nCount_RNA|nFeature_RNA|time|celltype|technology|dataset|predicted.id"))
  metadata$celltype <- NULL
  colnames(metadata)[colnames(metadata) == "predicted.id"] <- "celltype"
  seurat_obj@meta.data <- metadata
}

############################################################
# *. INFORMATION ABOUT THE INTEGRATED DATASET
############################################################

descriptive <- function(seurat_obj) {
  table(seurat_obj@celltype)
  table(seurat_obj@time)
  table(seurat_obj@dataset)
}

############################################################
# *. EVALUATION OF THE INTEGRATION
############################################################

# RLE PLOTS

plot_rle <- function(seurat_obj, 
                     assay = "integrated",
                     layer = "data",
                     n_genes = 5000,
                     n_cells = 150,
                     by_dataset = FALSE,
                     y_limits = if (by_dataset) c(-0.2, 0.2) else c(-0.65, 0.65),
                     seed = if (by_dataset) 34 else 45,
                     output_path = NULL) {

  set.seed(seed)
  
  gene_names <- rownames(seurat_obj[[assay]])
  cell_names <- colnames(seurat_obj[[assay]])

  sample_genes <- if (length(gene_names) <= n_genes) gene_names else sample(gene_names, n_genes)
  sample_cells <- if (length(cell_names) <= n_cells) cell_names else sample(cell_names, n_cells)
  
  data <- as.matrix(GetAssayData(seurat_obj, assay = assay, layer = layer))[sample_genes, sample_cells]

  medians <- apply(data, 1, median, na.rm = TRUE)
  rle_values <- sweep(log1p(data), 1, medians, "-")
  rle_df <- reshape2::melt(rle_values)
  colnames(rle_df) <- c("Gene", "Cell", "RLE")
  rle_df <- rle_df[is.finite(rle_df$RLE), ]
  
  if (by_dataset) {
    sample_metadata <- seurat_obj@meta.data[sample_cells, ]
    rle_df$Dataset <- sample_metadata[rle_df$Cell, "dataset"]
  }
  
  p <- ggplot(rle_df, aes_string(x = if (by_dataset) "Dataset" else "Cell", 
                                 y = "RLE",
                                 fill = if (by_dataset) "Dataset")) +
    geom_boxplot(outlier.shape = NA, show.legend = !by_dataset) +
    labs(title = paste("RLE plot", if (by_dataset) "by Dataset" else "after Integration"),
         y = "Relative Log Expression",
         x = NULL) +
         ylim(y_limits) +
         theme_minimal() +
         theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  if (by_dataset) {
    p <- p + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
      scale_fill_brewer(palette = "Set3")
  } else {
    p <- p + 
      theme(axis.text.x = element_blank())
  }
  if (!is.null(output_path)) {
    ggsave(output_path, plot = p, dpi = 1080, height = 7, width = if (by_dataset) 8 else 6)
  }
  return(p)
}

# UMAP PLOTS

plot_integration_umap <- function(seurat_obj,
                                  dims = 1:20,
                                  group_var = "celltype",
                                  split_var = "dataset",
                                  label = TRUE,
                                  label_size = 3,
                                  output_dir = "~/R/Figures/",
                                  prefix = "integrated") {
  
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj, verbose = FALSE, npcs = max(dims))
  seurat_obj <- RunUMAP(seurat_obj, dims = dims)
  
  plots <- list()
  
  plots$split_umap <- DimPlot(seurat_obj,
                              group.by = group_var,
                              split.by = split_var,
                              label = label,
                              label.size = label_size) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 9)) +
    labs(x = "UMAP 1", y = "UMAP 2")
  
  plots$combined_umap <- DimPlot(seurat_obj, group.by = group_var) + ggtitle("Celltypes")
  
  ggsave(file.path(output_dir, paste0(prefix, "_umap_split.png")),
         plot = plots$split_umap,
         dpi = 1080, height = 7, width = 12)
  
  ggsave(file.path(output_dir, paste0(prefix, "_umap_combined.png")),
         plot = plots$combined_umap,
         dpi = 1080, height = 7, width = 10)
  
  return(plots)
}

evaluate_integration_lisi <- function(seurat_obj,
                                 batch_var = "dataset",
                                 celltype_var = "celltype",
                                 time_var = "time",
                                 reduction = "umap",
                                 dims = 1:20,
                                 perplexity = 30,
                                 lisi_threshold = NULL,
                                 run_preprocessing = TRUE,
                                 plot = TRUE,
                                 verbose = TRUE) {
  
  # Check and run preprocessing if needed
  if (run_preprocessing) {
    if (verbose) message("Checking preprocessing requirements...")
    
    if (is.null(GetAssayData(seurat_obj, slot = "scale.data"))) {
      if (verbose) message("Running ScaleData...")
      seurat_obj <- ScaleData(seurat_obj)
    }
    
    if (!"pca" %in% names(seurat_obj@reductions)) {
      if (verbose) message("Running PCA...")
      seurat_obj <- RunPCA(seurat_obj, npcs = max(dims))
    }
    
    if (!reduction %in% names(seurat_obj@reductions)) {
      if (verbose) message("Running ", toupper(reduction), "...")
      seurat_obj <- RunUMAP(seurat_obj, dims = dims)
    }
  }
  
  # Calculate LISI scores
  embeddings <- Embeddings(seurat_obj, reduction = reduction)
  metadata <- seurat_obj@meta.data
  
  lisi_scores <- lisi::compute_lisi(
    X = embeddings,
    meta_data = data.frame(
      batch = metadata[[batch_var]],
      celltype = metadata[[celltype_var]]
    ),
    label_colnames = c("batch", "celltype"),
    perplexity = perplexity
  )
  
  lisi_df <- data.frame(
    Batch_LISI = lisi_scores$batch,
    CellType_LISI = lisi_scores$celltype
  )
  
  # Calculate summary statistics
  summary_stats <- data.frame(
    Metric = c("Batch LISI", "CellType LISI"),
    Mean = c(mean(lisi_df$Batch_LISI), mean(lisi_df$CellType_LISI)),
    Median = c(median(lisi_df$Batch_LISI), median(lisi_df$CellType_LISI)),
    SD = c(sd(lisi_df$Batch_LISI), sd(lisi_df$CellType_LISI)),
    N_Cells = nrow(lisi_df)
  )
  
  # Check if batch adjustment is needed
  n_datasets <- length(unique(metadata[[batch_var]]))
  threshold <- if (is.null(lisi_threshold)) n_datasets/2 else lisi_threshold
  
  if (mean(lisi_df$Batch_LISI) < threshold) {
    message("\n Warning! The LISI batch (", round(mean(lisi_df$Batch_LISI), 2), 
            ") doesn't indicate a good integration (", threshold, ").\n",
            "Please consider this steps in the downstream analysis:\n",
            "1. Adjust for cell proportions\n",
            "2. Adjust for technical parameters")
    
    # Calculate cell proportions
    barcodes <- rownames(metadata)
    cell_proportions <- metadata %>%
      group_by(.data[[batch_var]], .data[[time_var]], .data[[celltype_var]]) %>%
      summarise(count = n(), .groups = "drop") %>%
      group_by(.data[[batch_var]], .data[[time_var]]) %>%
      mutate(prop = count / sum(count)) %>%
      select(.data[[batch_var]], .data[[time_var]], .data[[celltype_var]], prop)
    
    # Add proportions to metadata
    seurat_obj@meta.data <- metadata %>%
      left_join(cell_proportions, 
                by = c(batch_var, time_var, celltype_var))
    
    rownames(seurat_obj@meta.data) <- barcodes
    colnames(seurat_obj) <- barcodes
  }
  
  # Generate plots if requested
  plots <- list()
  if (plot) {
    plots$batch <- ggplot(lisi_df, aes(x = Batch_LISI)) +
      geom_density(fill = "salmon", alpha = 0.7) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
      labs(title = paste("LISI - Batch (", batch_var, ")"),
           x = "LISI Score", y = "Density") +
      theme_minimal() +
      xlim(c(1, max(lisi_df$Batch_LISI)))
    
    plots$celltype <- ggplot(lisi_df, aes(x = CellType_LISI)) +
      geom_density(fill = "steelblue", alpha = 0.7) +
      labs(title = paste("LISI - Cell Type (", celltype_var, ")"),
           x = "LISI Score", y = "Density") +
      theme_minimal() +
      xlim(c(1, max(lisi_df$CellType_LISI)))
    
    plots$combined <- plots$batch + plots$celltype +
      plot_annotation(
        title = "Integration Quality Evaluation",
        subtitle = paste("Number of cells:", nrow(embeddings))
      )
    
    print(plots$combined)
  }
  
  # Return results
  return(list(
    seurat_obj = seurat_obj,
    lisi_scores = lisi_df,
    summary_stats = summary_stats,
    plots = if (plot) plots else NULL,
    threshold_check = mean(lisi_df$Batch_LISI) < threshold
  ))
}

############################################################
# *. CELL CYCLE SCORE
############################################################

plot_cell_cycle <- function(seurat_obj,
                            reduction = "pca",
                            group_var = "Phase",
                            split_var = "time",
                            dims = 1:20,
                            run_pca = FALSE,
                            output_path = "~/R/Figures/CC_plot.png",
                            plot_height = 7,
                            plot_width = 9,
                            dpi = 1080) {

  if (!"RNA" %in% names(seurat_obj@assays) || 
    is.null(GetAssayData(seurat_obj, slot = "scale.data"))) {
    message("Running ScaleData...")
    seurat_obj <- ScaleData(seurat_obj)
  }

  if (run_pca || !"pca" %in% names(seurat_obj@reductions)) {
    message("Running PCA...")
    seurat_obj <- RunPCA(seurat_obj, npcs = max(dims), verbose = FALSE)
  }

  cc_plot <- DimPlot(seurat_obj,
                     reduction = reduction,
                     group.by = group_var,
                     split.by = split_var) +
                     theme_minimal() +
                     labs(title = "Cell Cycle Phase Distribution",
                     subtitle = paste("Split by", split_var))

  if (!is.null(output_path)) {
    dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)
    ggsave(output_path, 
           plot = cc_plot,
           height = plot_height,
           width = plot_width,
           dpi = dpi)
    message("Plot saved to: ", output_path)
  }
  
  return(cc_plot)
}

############################################################
# *. FILTER ECM GENES
############################################################

filter_ecm_genes <- function(seurat_obj, ecm_genes_path) {
  ECM <- read_excel(ecm_genes_path)
  genes_to_keep <- ECM$Gene
  rownames_seurat <- rownames(seurat_obj)
  rows_to_keep <- which(rownames_obj %in% genes_to_keep)
  seurat_integrated_ME <- seurat_obj[rows_to_keep, ]
  paste("ECM Genes present in the dataset:", nrow(seurat_integrated_ME))
  return(seurat_integrated_ME)
} 

############################################################
# 5. DIFFERENTIAL EXPRESSION ANALYSIS
############################################################

DE_Analysis <- function(seurat_obj, 
                        condition, 
                        control, 
                        cell.prop = NULL,
                        logfc_thresh = 0, 
                        min_pct = 0.05,
                        output_file = NULL,
                        plot_path = NULL, 
                        plot_title = NULL) {
  
  Idents(seurat_obj) <- "time"
  
  # Set parameters
  params <- list(
    object = seurat_obj,
    ident.1 = condition,
    ident.2 = control,
    logfc.threshold = logfc_thresh,
    min.pct = min_pct,
    test.use = "MAST",
    layer = "data"
  )
  
  if (!is.null(cell.prop)) {
    params$latent.vars <- cell.prop
  }
  
  # Differential expression analysis
  DEG <- do.call(FindMarkers, params)
  DEG$Gene <- rownames(DEG)
  DEG$p_val_adj <- ifelse(DEG$p_val_adj == 0, 1e-305, DEG$p_val_adj)
  
  # Title
  if (is.null(plot_title)) {
    plot_title <- paste("Differential expression:", condition, "vs", control)
  }
  
  volcanoplot <- EnhancedVolcano(
    DEG,
    lab = DEG$Gene,
    x = "avg_log2FC",
    y = "p_val_adj",
    pCutoff = 0.05,
    FCcutoff = 0.75,
    title = plot_title,
    drawConnectors = TRUE,
    labSize = 3
  )
  
  if (!is.null(output_file)) {
    write.xlsx(DEG, output_file)
  }
  
  if (!is.null(plot_path)) {
    ggsave(plot_path, 
           plot = volcanoplot, 
           dpi = 1080, 
           height = 7, 
           width = 9)
  }
  
  # Results
  return(list(
    DEG_results = DEG,
    volcano_plot = volcanoplot
  ))
}

############################################################
# 6. PATHWAY ENRICHMENT ANALYSIS
############################################################

run_enrichment <- function(DEG_results,
                           pval_threshold = 0.05,
                           lfc_threshold = 0.75,
                           genes_already_filtered = FALSE,
                           organism_db = org.Mm.eg.db,
                           ontology = "BP",
                           keywords = c("neuro", "inflammation", "Wnt", "Rho", 
                                        "regression", "remodelling", "SCI",
                                        "axonal", "axon", "BMP", "spinal"),
                           exclude_terms = "muscle",
                           output_file = NULL) {
  
  # Verify required columns exist
  if (!all(c("p_val_adj", "avg_log2FC") %in% colnames(DEG_results))) {
    stop("The data frame must contain: 'p_val_adj' y 'avg_log2FC' columns")
  }
  
  # Filter genes if not already filtered
  if (!genes_already_filtered) {
    sig_genes <- DEG_results %>%
      filter(p_val_adj < pval_threshold & abs(avg_log2FC) > lfc_threshold)
    
    if (nrow(sig_genes) == 0) {
      stop("There are no significant genes with that thresholds: (p_val_adj < ", 
           pval_threshold, " y |log2FC| > ", lfc_threshold, ")")
    }
    gene_symbols <- rownames(sig_genes)
  } else {
    message("Genes already filtered")
    gene_symbols <- rownames(DEG_results)
  }
  
  # Convert gene symbols to ENTREZID
  gene_convert <- clusterProfiler::bitr(gene_symbols,
                                        fromType = "SYMBOL",
                                        toType = "ENTREZID",
                                        OrgDb = organism_db)
  
  if (nrow(gene_convert) == 0) {
    stop("Error in SYMBOL to ENTREZID")
  }
  
  # GO enrichment analysis
  go_enrich <- enrichGO(
    gene = gene_convert$ENTREZID,
    OrgDb = organism_db,
    keyType = "ENTREZID",
    ont = ontology,
    readable = TRUE
  )
  
  # Filter terms of interest
  if (!is.null(keywords) && nrow(go_enrich) > 0) {
    keyword_pattern <- paste(keywords, collapse = "|")
    matched_terms <- go_enrich@result %>%
      filter(str_detect(Description, regex(keyword_pattern, ignore_case = TRUE)))
    
    if (!is.null(exclude_terms)) {
      exclude_pattern <- paste(exclude_terms, collapse = "|")
      matched_terms <- matched_terms %>%
        filter(!str_detect(Description, regex(exclude_pattern, ignore_case = TRUE)))
    }
    
    if (nrow(matched_terms) > 0) {
      go_enrich@result <- matched_terms
    } else {
      warning("0 GO terms")
      return(NULL)
    }
  }
  
  # Generate plot
  if (nrow(go_enrich@result) > 0) {
    enrich_plot <- barplot(go_enrich, showCategory = 15) +
      theme_minimal() +
      theme(
        axis.text.y = element_text(size = 9),
        plot.caption = element_text(hjust = 0)
      ) +
      labs(
        title = "GO Pathway Enrichment",
        caption = if (!is.null(keywords)) {
          paste("Terms related with:", paste(keywords, collapse = ", "))
        } else {
          "GO significative terms"
        }
      )
    
    # Save plot if output path is provided
    if (!is.null(output_file)) {
      ggsave(output_file, 
             plot = enrich_plot,
             dpi = 1080,
             height = 7,
             width = 9)
    }
    
    return(list(
      enrich_result = go_enrich,
      plot = enrich_plot,
      gene_conversion = gene_convert
    ))
  } else {
    warning("No results to show")
    return(NULL)
  }
}

############################################################
# 7. SPECIFIC CELL TYPE DIFFERENTIAL EXPRESSION
############################################################

HeatMap <- function(seurat_obj, 
                    ident_column = "time",
                    condition = "7dpi",
                    control = "Uninjured",
                    celltypes = c("Astrocytes", "Neurons", "Microglia", "OPCs", "Oligodendrocytes", "Endothelial"),
                    genes_up = NULL,
                    genes_down = NULL,
                    logfc_threshold = 0,
                    min_pct = 0.05,
                    test_use = "MAST",
                    output_file = "~/R/Figures/Heatmap_Acute.png",
                    plot_width = 14,
                    plot_height = 7,
                    plot_dpi = 1080) {
  
  # Set identity classes
  Idents(seurat_obj) <- ident_column
  
  # Initialize empty list to store DEG results
  deg_results <- list()
  
  # Find DEGs for each cell type
  for (ct in celltypes) {
    # Subset the Seurat object
    cells <- subset(seurat_obj, subset = celltype == ct)
    
    # Find markers
    deg <- FindMarkers(cells, 
                       ident.1 = condition, 
                       ident.2 = control, 
                       logfc.threshold = logfc_threshold, 
                       min.pct = min_pct, 
                       test.use = test_use)
    
    # Add gene and celltype information
    deg$Gene <- rownames(deg)
    deg$celltype <- ct
    
    # Store results
    deg_results[[ct]] <- deg
  }
  
  # Combine all results
  combined_results <- bind_rows(deg_results)
  
  # Filter genes if provided
  if (!is.null(genes_up) | !is.null(genes_down)) {
    combined_results <- combined_results %>%
      filter(Gene %in% genes_down | Gene %in% genes_up)
  }
  
  # Prepare heatmap data
  heatmap_data <- combined_results %>%
    select(Gene, celltype, avg_log2FC) %>%
    pivot_wider(names_from = celltype, values_from = avg_log2FC) %>%
    column_to_rownames("Gene") %>%
    as.matrix()
  
  # Replace NA with 0
  heatmap_data[is.na(heatmap_data)] <- 0
  
  # Create heatmap
  heatmap <- Heatmap(
    heatmap_data,
    col = circlize::colorRamp2(c(-7, 0, 7), c("blue", "white", "red")),
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 8),
    name = "log2FC"
  )
  
  # Save plot
  if (!is.null(output_file)) {
    ggsave(output_file, heatmap, dpi = plot_dpi, height = plot_height, width = plot_width)
  }
  
  # Return the heatmap
  return(heatmap)
}


############################################################
# APPLICATION OF THE PIPELINE
############################################################


