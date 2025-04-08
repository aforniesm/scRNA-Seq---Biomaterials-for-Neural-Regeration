############################################################
# 0. R PACKAGES NEEDED
############################################################

if (!require("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager")

# If the packages are not installed use this code:
# BiocManager::install(c("glmGamPoi", "biomaRt", "EnhancedVolcano", "speckle", "lisi", "AnnotationDbi","org.Mm.eg.db"))
# install.packages(c("Seurat", "Matrix", "tidyverse", "ggvenn", "tibble", "circlize", "patchwork))

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
library(patchwork)
library(readxl)
library(openxlsx)
library(clusterProfiler)
library(enrichplot)

############################################################
# 1. LOAD THE DATA 
############################################################

# As the datasets of single cell RNA-Seq have multiple formats, and each study employs different
# nomenclatures, you should load your datasets as Seurat objects. Thank you!

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
  return(seurat_obj)
}

# Create a list with the Seurat Objects that you want to integrate
# seurat_list <- list()
  
############################################################
# 4. INTEGRATION
############################################################  

prepare_integration <- function(seurat_list) {
  seurat_list <- lapply(X = seurat_list, FUN = function(x) {x <- ScaleData(x, features =  features)
  x <- RunPCA(x, features = features, npcs = 30)
  x <- RunUMAP(x, dims = 1:20)})
  return(seurat_list)
  }

integration <- function(seurat_list, features) {
  anchors_integration <- FindIntegrationAnchors(object.list = seurat_list, dims = 1:20, verbose = TRUE, reduction = "rpca", anchor.features = features, scale = F, normalization.method = "LogNormalize")
  seurat_integrated <- IntegrateData(anchorset = anchors_integration, dims = 1:20, verbose = TRUE, normalization.method = "LogNormalize")
  return(seurat_integrated)
}

############################################################
# *. UMAP BEFORE INTEGRATION (Add more maps if you integer more datasets)
############################################################  

umap_visualization <- function(seurat_list, path) {
  
  plot_list <- list()
  for (i in seq_along(seurat_list)) {
    plot_list[[i]] <- DimPlot(
      seurat_list[[i]],
      split.by = "dataset",
      group.by = "predicted.id",
      label = TRUE,
      label.size = 3
    ) +
      theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7)) +
      labs(x = "UMAP 1", y = "UMAP 2")
  }
  
  combined_plot <- Reduce(`+`, plot_list)
  
  ggsave(path, combined_plot, dpi = 1080, height = 14, width = 8)
}

############################################################
# *. METADATA
############################################################  

metadata_cleaning <- function(seurat_obj) {
  metadata <- seurat_obj@meta.data
  selected_columns <- c("nCount_RNA", "nFeature_RNA", "time", "celltype", "technology", "dataset", "predicted.id")
  metadata <- metadata[, selected_columns, drop = FALSE]
  metadata$celltype <- NULL
  colnames(metadata)[colnames(metadata) == "predicted.id"] <- "celltype"
  seurat_obj@meta.data <- metadata
  return(seurat_obj)
}

metadata_cleaning2 <- function(seurat_obj, reference_dataset) {
  metadata <- seurat_obj@meta.data
  selected_columns <- c("nCount_RNA", "nFeature_RNA", "time", "celltype", "technology", "dataset", "predicted.id")
  existing_columns <- selected_columns[selected_columns %in% colnames(metadata)]
  metadata <- metadata[, existing_columns, drop = FALSE]
  metadata$celltype <- ifelse(metadata$dataset == reference_dataset & is.na(metadata$predicted.id),
                              metadata$celltype,  
                              metadata$predicted.id)  
  metadata$predicted.id <- NULL
  seurat_obj@meta.data <- metadata
  return(seurat_obj)
}

############################################################
# *. INFORMATION ABOUT THE INTEGRATED DATASET
############################################################

descriptive <- function(seurat_obj) {
  df_celltype <- as.data.frame(table(seurat_obj@meta.data$celltype))
  df_time <- as.data.frame(table(seurat_obj@meta.data$time))
  df_dataset <- as.data.frame(table(seurat_obj@meta.data$dataset))
  return(list(celltype = df_celltype, time = df_time, dataset = df_dataset))
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
# *. CELL PROPORTION
############################################################

celltype_proportions <- function(seurat_obj, dataset_var = "dataset", time_var = "time", celltype_var = "celltype") {
  metadata <- seurat_obj@meta.data
  barcodes <- rownames(metadata)
  
  cell_proportions <- metadata %>%
    dplyr::group_by(.data[[dataset_var]], .data[[time_var]], .data[[celltype_var]]) %>%
    dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
    dplyr::group_by(.data[[dataset_var]], .data[[time_var]]) %>%
    dplyr::mutate(prop = count / sum(count)) %>%
    dplyr::select(
      !!dataset_var := .data[[dataset_var]],
      !!time_var := .data[[time_var]],
      !!celltype_var := .data[[celltype_var]],
      prop
    )
  
  prop_table <- cell_proportions %>%
    tidyr::pivot_wider(
      names_from = .data[[celltype_var]],
      values_from = prop,
      values_fill = 0
    )
  
  seurat_obj@meta.data <- metadata %>%
    dplyr::left_join(
      prop_table,
      by = setNames(c(dataset_var, time_var), c(dataset_var, time_var))
    )
  
  rownames(seurat_obj@meta.data) <- barcodes
  return(seurat_obj)
}


############################################################
# *. CELL CYCLE SCORE
############################################################

plot_cell_cycle <- function(seurat_obj, 
                            cell_cycle_markers_path = "~/R/Data/Mus_musculus.csv",
                            output_path = NULL,
                            reduction = "pca",
                            group_var = "Phase",
                            split_var = "time",
                            plot_height = 7,
                            plot_width = 9,
                            dpi = 1080) {
  
  cell_cycle_markers <- read.csv(cell_cycle_markers_path)
  
  s_genes <- cell_cycle_markers %>% dplyr::filter(phase == "S") %>% pull("geneID")
  g2m_genes <- cell_cycle_markers %>% dplyr::filter(phase == "G2/M") %>% pull("geneID")
  
  s_genes <- AnnotationDbi::select(org.Mm.eg.db, 
                                   keys = s_genes, 
                                   keytype = "ENSEMBL", 
                                   columns = c("SYMBOL"))
  g2m_genes <- AnnotationDbi::select(org.Mm.eg.db, 
                                     keys = g2m_genes, 
                                     keytype = "ENSEMBL", 
                                     columns = c("SYMBOL"))
  
  seurat_obj <- CellCycleScoring(seurat_obj, 
                                 g2m.features = g2m_genes$SYMBOL, 
                                 s.features = s_genes$SYMBOL)
  
  cc_score <- DimPlot(seurat_obj, 
                      reduction = reduction, 
                      group.by = group_var, 
                      split.by = split_var)
  
  ggsave(output_path, cc_score, dpi = dpi, height = plot_height, width = plot_width)
  
  message("Plot saved to: ", output_path)
  return(seurat_obj)
}


############################################################
# *. FILTER ECM GENES
############################################################

filter_ecm_genes <- function(seurat_obj, ecm_genes_path) {
  ECM <- read_excel(ecm_genes_path)
  genes_to_keep <- ECM$Gene
  rownames_seurat <- rownames(seurat_obj)
  rows_to_keep <- which(rownames_seurat %in% genes_to_keep)  # Cambio aquí
  seurat_integrated_ME <- seurat_obj[rows_to_keep, ]
  message("ECM Genes present in the dataset:", nrow(seurat_integrated_ME))
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
    cell.prop <- c("Astrocytes", "Neurons", "OPCs", "Oligodendrocytes", "Microglia", "Endothelial")
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


#####################################################################################################
################################ APPLICATION OF THE PIPELINE ########################################
#####################################################################################################

# 1. Load the data as Seurat Objects
load("~/R/Data/Seurats_For_Acute.RData")

# 2. Normalization
seurat_norm <- normalize_data(seurat_obj = seurat)
seurat2_norm <- normalize_data(seurat_obj = seurat2)
seurat3_norm <- normalize_data(seurat_obj = seurat3)
seurat5_norm <- normalize_data(seurat_obj = seurat5)
rm(seurat, seurat2, seurat3, seurat5)

# 3. Cell Annotation
seurat_annot <- cell_annotation(seurat_obj = seurat_norm, reference = seurat5_norm)
seurat2_annot <- cell_annotation(seurat_obj = seurat2_norm, reference = seurat5_norm)
seurat3_annot <- cell_annotation(seurat_obj = seurat3_norm, reference = seurat5_norm)
rm(seurat_norm, seurat2_norm, seurat3_norm, seurat5_norm)

# 4. Integration
seurat_list <- list(seurat_annot, seurat2_annot, seurat3_annot)
rm(seurat_annot, seurat2_annot, seurat3_annot)
features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 3500)
seurat_list_prep <- prepare_integration(seurat_list)
# features <- c(features, genes_of_interest) For sub-chronic or specific gene
seurat_integrated <- integration(seurat_list = seurat_list_prep, features = features)
Idents(seurat_integrated) <- "time"

# 5. UMAP before Integration
umap_before_integration <- umap_visualization(seurat_list = seurat_list, path = "~/R/Figures/UMAP_Before_Integration.png")

# 6. Post-Integration Steps
# Metadata
seurat_integrated <- metadata_cleaning(seurat_integrated)
# Descriptive
descriptive(seurat_obj = seurat_integrated)
# Preparation of evaluation of integration
seurat_integrated <- ScaleData(seurat_integrated)
seurat_integrated <- RunPCA(seurat_integrated, verbose = FALSE)
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:20)
# UMAP
plot_integration_umap(seurat_obj = seurat_integrated, prefix = "Acute")
# RLE
plot_rle(seurat_obj = seurat_integrated, output_path = "~/R/Data/RLE-Acute.png")
plot_rle(seurat_obj = seurat_integrated, by_dataset = T, output_path = "~/R/Data/RLE-Acute(Dataset).png")
# LISI
evaluate_integration_lisi(seurat_obj = seurat_integrated, lisi_threshold = 2)
# Cell proportion calculation
seurat_integrated <- celltype_proportions(seurat_obj = seurat_integrated)

# 7. Cell cycle 
seurat_integrated <- plot_cell_cycle(seurat_obj = seurat_integrated, output_path = "~/R/Figures/CC_Acute.png")

# 8. ECM genes
seurat_ME <- filter_ecm_genes(seurat_obj = seurat_integrated, ecm_genes_path = "~/R/Data/MAIN_ECM_Proteome.xlsx")
seurat_ME@meta.data <- seurat_ME@meta.data %>%
  rename(`Endothelial cells` = "Endothelial")

# 9. DGE
DEG <- DE_Analysis(seurat_obj = seurat_ME, condition = "7dpi", control = "Uninjured", cell.prop = T, plot_path = "~/R/Figures/Volcano_Acute.png", plot_title = "Differential expression analysis in Acute vs Uninjured conditions", output_file = "~/R/Data/DEG_Acute_7dpi.xlsx" )

# 10. Enrichment analysis
run_enrichment(DEG_results = DEG, output_file = "~/R/Figures/Go_Enrich_Acute.png")

# 11. Cell type DGE
HeatMap(seurat_obj = seurat_ME)
