## This code contains a series of functions that can be used to process AAV SPLiTseq data
##   that has gone through the STEP1 and STEP2 AAV_SPLiTseq R scripts.

## Analysis of Single cell RNAseq data analysis requires significant customization based on the
##    data being analyzed and the questions being asked. For that reason, the following code
##    provides a series of functions that can be used as tools to perform some common
##    Seurat based single cell analyses and data manipulations. These should be considered
##    a possible starting point for such analyses. That said, anyone familiar with
##    single cell RNAseq analysis may find that this is the point in the pipeline where
##    it is more profitable to rely on their own expertise.


suppressPackageStartupMessages({
  library(glmGamPoi)
  library(harmony)
  library(readr)
  library(cowplot)
  library(details)
  library(viridis)
  library(scCustomize)
  library(patchwork)
  library(Seurat)
  library(tidyseurat)
  library(dplyr)
  library(tidyverse)
  
  # Doublet finder dependencies
  library(Matrix)
  library(fields)
  library(KernSmooth)
  library(ROCR)
  library(parallel)
  library(remotes)
  ## Use the lines below to install DoubletFinder if needed.
  ##    EDIT BCL: I had to manually change the file download method. Otherwise, I was getting
  ##              a 'cannot open URL' message. The options() code below fixed it.
  # options(download.file.method = "wininet")
  # remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
  library(DoubletFinder)
  
  ## Seurat v5 documentation says they use presto to dramatically improve the speed of DE cluster
  ##    analysis. You can install it using the code below:
  # install.packages("devtools")
  #devtools::install_github("immunogenomics/presto")
  library(presto)
  
  ## Other
  library(officer)
  library(rvg)
  library(writexl)
  library(ggplotify)
  
})








##  step3_seurat_preprocessing(seurat_data)
##
##  Run this function ONCE on your seurat object from STEP2. This just creates a new column
##      in the meta data that will store notes on all modifications made throughout the
##      STEP3 script.
##  Important Note: If your seurat_data already has a Notes column in it's metadata, this function
##      will over-write it. So, make sure you don't run this on a seurat object that you have already
##      started subsetting using the STEP3 script.

step3_seurat_preprocessing <- function(seurat_data) {
  DefaultAssay(seurat_data) <- 'RNA'
  seurat_data@meta.data$Notes <- NA
  seurat_data@meta.data$Notes[1] <- 'This column tracks edits made during cell annotation'
  return(seurat_data)
}













##  Make_QC_vs_Cluster_Plots(seurat_data, features_to_plot = c("nFeature_RNA", "nCount_RNA", "percent_mito"), 
##                           group_by_id="seurat_clusters")
##    features_to_plot is optional. It will default to c("nFeature_RNA", "nCount_RNA", "percent_mito").
##    group_by_id is optional. It will default to 'seurat_clusters'.

Make_QC_vs_Cluster_Plots <- function(seurat_data, features_to_plot = c("nFeature_RNA", "nCount_RNA", "percent_mito"), 
                                     group_by_id="seurat_clusters") {
  output_plot <- VlnPlot(seurat_data, group.by=group_by_id, 
          features = features_to_plot, 
          ncol = 1, pt.size=0, log = FALSE) 
  
  return(output_plot)
  output_plot
}




## pull_cluster_markers(seurat_data, dataset_prefix=data_prefix, save_directory=save_dir, 
##                      cluster_id='seurat_clusters', save_data_to_csv = TRUE,
##                      use_default_options=TRUE, custom_options=c())
## 
## The default options are:
##          marker_test <- 'wilcox'
##          pos_neg_option <- TRUE
##          min_pct_val <- 0.01
##          logfc_threshold_val <- 0.001
## To use custom options, provide an array with the values you would like to use in the same order
##    as the default options listed above. e.g., here's how you would feed the default values above
##    into the code as custom options:  custom_options = c('wilcox',TRUE,0.01,0.001)

pull_cluster_markers <- function(seurat_data, data_prefix, save_directory=save_dir, 
                                 cluster_id='seurat_clusters', save_data_to_csv = TRUE,
                                 use_default_options=TRUE, custom_options=c()) {
  
  if (use_default_options) {
    marker_test <- 'wilcox'
    pos_neg_option <- TRUE
    min_pct_val <- 0.01
    logfc_threshold_val <- 0.001
  } else {
    marker_test <- custom_options[1]
    pos_neg_option <- custom_options[2]
    min_pct_val <- custom_options[3]
    logfc_threshold_val <- custom_options[4]
  }
  
  # Find all markers
  options(future.globals.maxSize= 8100000000)
  
  Idents(seurat_data) <- cluster_id
  all_markers <- FindAllMarkers(
    seurat_data,
    test.use = marker_test,
    only.pos = pos_neg_option,
    min.pct = min_pct_val,
    logfc.threshold = logfc_threshold_val
  )
  
  all_markers <- all_markers %>%
    group_by(cluster) %>%
    arrange(desc(avg_log2FC)) %>%
    arrange(cluster)
  
  top10_markers <- all_markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)
  
  if (save_data_to_csv) {
    write_csv(all_markers, file = paste0(save_dir,'/',data_prefix,"_allclustermarkers.csv"))
    write_csv(top10_markers, file = paste0(save_dir,'/',data_prefix,"_top10markers.csv"))
  }
  
  top_markers_list <- list()
  top_markers_list[['all_cluster_markers']] <- all_markers
  top_markers_list[['top10_cluster_markers']] <- top10_markers
  return(top_markers_list)
  
}




## plot_DoHeatmap(seurat_data, cell_markers, group_by_id = 'seurat_clusters')
##
## This is basically just the DoHeatmap function. The only purpose of encapsulating it in
##    this function is to keep the format for calling it similar to that used for the other
##    plotting functions in this script.


plot_DoHeatmap <- function(seurat_data, cell_marker_list, group_by_id = 'seurat_clusters') {
  
  output_plot <- DoHeatmap(seurat_data, group.by ="seurat_clusters", features=as.character(unlist(cell_markers)))
  
  output_plot
  return(output_plot)
  
}








## make_celltype_dotplot(seurat_data, cell_marker_list)
##
##  This function takes a list of cell type markers with the cell type name
##      as the name of the list entry, and the entries themselves as vectors
##      of gene names.
##  Note that this plot will look best if you use short names for the cell types
##      long names are likely to overlap. This can be fixed by extending the plot
##      horizontally.
##
##  broad_cell_markers <- list(
##    'Cell_TypeA' = c('Gene1',	'Gene2',	'Gene3'),
##    'Cell_TypeB' = c('Gene4','Gene5'),
##    ...
##  )

make_celltype_dotplot <- function(seurat_data, cell_marker_list, group_by_id = 'seurat_clusters') {
  
  dot_plots <- list()
  
  cell_markers <- cell_marker_list[[1]]
  
  dot_plot <- DotPlot(object = seurat_data, group.by = group_by_id, features = cell_markers) + scale_x_discrete(guide = guide_axis(angle = 90)) +
    theme(legend.position = 'none', axis.title = element_blank()) +
    labs(title = names(cell_marker_list)[1]) + theme(plot.title = element_text(hjust = 0.5, size=12))
  
  dot_plots[[1]] <- dot_plot
  
  if (length(cell_marker_list) > 2) {
    
    for (i in 2:(length(cell_marker_list)-1)) {
      
      cell_markers <- cell_marker_list[[i]]
      
      dot_plot <- DotPlot(object = seurat_data, features = cell_markers) + scale_x_discrete(guide = guide_axis(angle = 90)) +
        theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
              axis.line.y = element_blank(), legend.position = 'none', axis.title = element_blank()) +
        labs(title = names(cell_marker_list)[i]) + theme(plot.title = element_text(hjust = 0.5, size=12))
      
      dot_plots[[i]] <- dot_plot
      
    }
    
    cell_markers <- cell_marker_list[[length(cell_marker_list)]]
    
    dot_plot <- DotPlot(object = seurat_data, features = cell_markers) + scale_x_discrete(guide = guide_axis(angle = 90)) +
      theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
            axis.line.y = element_blank(), axis.title = element_blank()) +
      labs(title = names(cell_marker_list)[length(cell_markers)]) + theme(plot.title = element_text(hjust = 0.5, size=12))
    
    dot_plots[[length(cell_marker_list)]] <- dot_plot
    
  } else {
    
    cell_markers <- cell_marker_list[[length(cell_marker_list)]]
    
    dot_plot <- DotPlot(object = seurat_data, features = cell_markers) + scale_x_discrete(guide = guide_axis(angle = 90)) +
      theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
            axis.line.y = element_blank(), axis.title = element_blank()) +
      labs(title = names(cell_marker_list)[length(cell_markers)]) + theme(plot.title = element_text(hjust = 0.5, size=12))
    
    dot_plots[[length(cell_marker_list)]] <- dot_plot
    
  }
  
  
  if (length(dot_plots) > 2) {
    
    dot_plot_patch <- dot_plots[[1]]+dot_plots[[2]]+plot_layout(ncol=2)
    for (i in 3:length(dot_plots)) {
      dot_plot_patch <- dot_plot_patch + dot_plots[[i]] + plot_layout(ncol=i)
    }
    
  } else {
    dot_plot_patch <- dot_plots[[1]]+dot_plots[[2]]+plot_layout(ncol=2)
  }
  
  dot_plot_patch
  return(dot_plot_patch)
}
























##  run_SingleR_cell_annotation(seurat_data, 
##                              reference_data, reference_data_notes,
##                              data_prefix, save_dir,
##                              convert_ensembl_to_hgnc = TRUE,
##                              remove_native_cells = FALSE,
##                              use_varFeatures_for_SingleR = TRUE,
##                              nFeatures_for_SingleR_comp = 3000,
##                              save_SingleR_predictions_dataframe = TRUE,
##                              save_SingleR_QC_figures = TRUE,
##                              save_seurat_SingleR = TRUE)
##
## This is a complex function. I've automated it as much as possible, but you'll need
##    to ensure that you provide it with the correct inputs:
##
## seurat_data - The seurat object you want annotated
## reference_data - The reference dataset you want to use. This MUST be a seurat object.
## reference_data_notes - There is no generic way for the script to know what reference data
##                        you are giving it. This variable should be a string containing whatever
##                        info on the reference dataset you want saved as part of the RunInfo output file.
##                        At a minimum, I'd suggest a citation for the article the data came from.
## data_prefix - The tag you want added to the saved filenames identifying the dataset. This should be short,
##               something like 'Retina'. The code will automatically add date stamps to the saved files
##               so they can be definitely identified as coming from the same run.
## save_dir - The directory where files will be saved. Do NOT add a '/' to the end of this variable.
## convert_ensembl_to_hgnc - Set to TRUE if the reference data has it's gene names in ENSEMBL format
## remove_native_cells - Set to TRUE if your reference data has a 'native cells' category and you want to remove it
## use_varFeatures_for_SingleR - Set to TRUE if you want to perform SingleR using only the top variable features
##                               from the experimental data. If FALSE, all genes will be used.
## nFeatures_for_SingleR_comp - Number of variable features to use if use_varFeatures_for_SingleR is set to TRUE
## save_SingleR_predictions_dataframe, save_SingleR_QC_figures, save_seurat_SingleR - Toggles to control the output
##      of the script. These toggles control the files that will be saved to your hard drive in the save_dir folder.
## 
## The function will return the seurat_data object with a SingleR_predictions column added to it's metadata.


run_SingleR_cell_annotation <- function(seurat_data, 
                                        reference_data, reference_data_notes,
                                        data_prefix, save_dir,
                                        convert_ensembl_to_hgnc = TRUE,
                                        remove_native_cells = FALSE,
                                        use_varFeatures_for_SingleR = TRUE,
                                        nFeatures_for_SingleR_comp = 3000,
                                        save_SingleR_predictions_dataframe = TRUE,
                                        save_SingleR_QC_figures = TRUE,
                                        save_seurat_SingleR = TRUE) {
  
  
  
  
  if (!use_varFeatures_for_SingleR) {
    nFeatures_for_SingleR_comp <- NA
  }
  
  
  ##  Optional removal of "native" cell category. It basically adds background noise to the cell type mapping
  ##    this would probably be useful if we were getting super high confidence mappings, but
  ##    that isn't the case with the Retina sample dataset. So, instead it is basically making
  ##    it more difficult to identify the higher confidence mappings because it reduces the
  ##    overall magnitude of the already relatively small values.
  
  if (remove_native_cells) {
    Idents(reference_data) <- 'cell_type'
    reference_data <- subset(reference_data, idents = c("native cell"), invert = TRUE)
  }
  
  
  
  
  ##################################################
  ##
  ##  Processing data from Reference Genome Sets to prepare it for use in SingleR
  ##
  ##  1) First, we need to create a SingleCellExperiment object for the reference data. This SCE
  ##        will be used for further processing of the reference set data and ultimately for
  ##        the SingleR analysis.
  ##  2) Filtering CellType_GeneRef_Set so it only includes genes present in the G_list
  ##        data frame of ensembl_gene_id to hgnc_symbol conversions
  ##  3) Pulling the ensembl_gene_ids from the filtered reference set then left_join'ing
  ##        them to the G_list to create a data frame with only the ensemble_gene_ids
  ##        and hgnc_symbols present in the filtered reference data and with those gene
  ##        names in the correct order
  ##  4) Replacing the ensembl gene names in the reference data with their hgnc_symbols
  
  ref_data_SCE <- as.SingleCellExperiment(reference_data, assay = 'RNA')
  
  
  if (convert_ensembl_to_hgnc) {
    
    ######################################
    ##
    ##  Pulling BioMart info to enable Ensembl <--> hgnc_symbol conversion
    
    mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
    avail_attributes <- listAttributes(mart)
    
    genes <- rownames(seurat_data[['RNA']])
    #G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
    G_list <- getBM(filters= "hgnc_symbol", attributes= c("ensembl_gene_id","hgnc_symbol"), values=genes, mart= mart)
    
    ## Filtering reference data so that it only includes entries found in the ensembl <--> hgnc conversion database
    ##    created based on the genes found in the experimental seurat data
    ref_data_SCE_Filt <- ref_data_SCE[rownames(ref_data_SCE) %in% G_list$ensembl_gene_id,]
    Ref_GeneIDs <- data.frame(ensembl_gene_id = rownames(ref_data_SCE_Filt))
    Ref_GeneIDs <- left_join(Ref_GeneIDs,G_list,by="ensembl_gene_id")
    rownames(ref_data_SCE_Filt) <- Ref_GeneIDs$hgnc_symbol
    
  } else if (convert_ensembl_to_hgnc == FALSE) {
    
    genes <- rownames(seurat_data[['RNA']])
    ref_data_SCE_Filt <- ref_data_SCE[rownames(ref_data_SCE) %in% genes,]
    
  }
  
  
  
  ##########################################
  ##
  ##  Running SingleR using all common genes in the experimental vs reference set data can result
  ##  in cell ID confidence scores that are relatively low. This is a variation to try
  ##    and improve the cell ID predictions.
  ##  Essentially, it runs SingleR using only the top variable features from the experimental dataset
  ##  
  
  if (use_varFeatures_for_SingleR) {
    
    seurat_data <- FindVariableFeatures(object = seurat_data, nfeatures = nFeatures_for_SingleR_comp)
    seurat_variable_features <- VariableFeatures(seurat_data)
    
    ## Filtering the processed reference data to match the list of variable features identified
    ##    in the experimental data
    ref_data_SCE_Filt <- ref_data_SCE_Filt[rownames(ref_data_SCE_Filt) %in% seurat_variable_features,]
    
  }
  
  
  ####################################################
  ##
  ##  Creating a SingleCellExperiment object for the experimental data
  ##
  
  my_SCE <- as.SingleCellExperiment(seurat_data, assay = 'RNA')
  
  ####################################################
  ##
  ##  Checking whether the gene IDs in the experimental SCE match the gene IDs in the
  ##    reference data set. This code is a check for the workflow if all experimental
  ##    genes are going to be used for SingleR.
  ##  This is currently commented out because it is not a relevant comparison for the
  ##    top variable features methodology.
  
  #my_SCE_GeneIDs <- data.frame(hgnc_symbol = rownames(my_SCE))
  #my_SCE_GeneIDs <- data.frame(hgnc_symbol = my_SCE_GeneIDs[as.character(my_SCE_GeneIDs$hgnc_symbol) %in% Retina_GeneIDs$hgnc_symbol,])
  
  #if (nrow(Retina_GeneIDs) == nrow(my_SCE_GeneIDs)) {
  #  print("Gene Lists Match")
  #}
  
  
  ####################################################
  ##
  ##  Filtering the experimental SCE so that it contains the same genes that are present
  ##    in the reference data. We already filtered the reference data so that it contains only
  ##    genes present in the experimental data, but it is possible there were genes in the
  ##    experimental data not present in the reference data. This additional filtering step
  ##    should ensure both the ref and exp data contain the same set of genes.
  ##
  
  my_SCE_Filt <- my_SCE[rownames(my_SCE) %in% rownames(ref_data_SCE_Filt),]
  
  if (nrow(my_SCE_Filt) != nrow(ref_data_SCE_Filt)) {
    print("Warning - Ref vs Exp data gene count mismatch")
  }
  
  
  
  ## Ref data needs to be log normalized prior to analysis with SingleR
  ## Experimental data does NOT need to be log normalized. Raw counts are fine for the experimental data.
  
  ref_data_SCE_Filt <- logNormCounts(ref_data_SCE_Filt)
  
  ## Now that we've got all the data properly filtered and processed, we can actually run the SingleR automated
  ##    cell annotation.
  SingleR_predictions <- SingleR(test=my_SCE_Filt, ref=ref_data_SCE_Filt, labels=ref_data_SCE_Filt$cell_type, assay.type.test = 1)
  
  
  if (save_SingleR_predictions_dataframe) {
    
    thetime <- Sys.time()
    Date_Tag <- paste0(substr(thetime,3,4),substr(thetime,6,7),substr(thetime,9,10),substr(thetime,8,8),
                       substr(thetime,12,13),'h',substr(thetime,15,16),'m')
    SingleR_output_name <- paste0("SingleR_Predictions-",data_prefix,"-",Date_Tag)
    
    ## Saving SingleR predictions data structure
    saveRDS(SingleR_predictions, paste0(save_dir,"/",SingleR_output_name,".rds"))
    
  }
  
  
  
  if (save_SingleR_QC_figures) {
    
    figures_PP <- read_pptx()
    windowsFonts("Times" = windowsFont("Times"))
    text_mod <- fp_text_lite(font.size=22, color='black')
    
    SingleR_table <- as.data.frame(table(SingleR_predictions$labels))
    SingleR_table$Pct <- SingleR_table$Freq / sum(SingleR_table$Freq) * 100
    SingleR_table$Pct <- round(SingleR_table$Pct,1)
    colnames(SingleR_table) <- c('Cell_Type','Count','Percent')
    ft_SingleR <- flextable(SingleR_table)
    ft_SingleR <- autofit(ft_SingleR)
    
    Refset_table <- as.data.frame(table(reference_data@meta.data$cell_type))
    Refset_table$Pct <- Refset_table$Freq / sum(Refset_table$Freq) * 100
    Refset_table$Pct <- round(Refset_table$Pct,1)
    colnames(Refset_table) <- c('Cell_Type','Count','Percent')
    ft_Refset <- flextable(Refset_table)
    ft_Refset <- autofit(ft_Refset)
    
    slide_title <- 'Reference Set Distribution (left)\nSingleR Predictions (right)'
    
    figures_PP <- add_slide(figures_PP,layout = "Title and Content", master = "Office Theme")
    figures_PP <- ph_with(figures_PP, value = slide_title, location = ph_location_type(type = "title"))
    figures_PP <- ph_with(figures_PP, value = ft_Refset, ph_location_left())
    figures_PP <- ph_with(figures_PP, value = ft_SingleR, ph_location_right())
    
    slide_title <- paste0('SingleR Prediction Score Heatmap')
    normal_plot <- as.ggplot(plotScoreHeatmap(SingleR_predictions))  
    #vectorized_plot <- rvg::dml(ggobj = normal_plot)
    figures_PP <- add_slide(figures_PP,layout = "Title and Content", master = "Office Theme")
    figures_PP <- ph_with(figures_PP, value = slide_title, location = ph_location_type(type = "title"))
    figures_PP <- ph_with(figures_PP, value = normal_plot, ph_location_type(type = "body"))
    
    slide_title <- paste0('SingleR Metric: delta median')
    normal_plot <- plotDeltaDistribution(SingleR_predictions, size=NA, show='delta.med') + theme(legend.position = 'none') + plot_annotation(title='SingleR Metric: delta median', theme = theme(plot.title = element_text(hjust = 0.5)))
    vectorized_plot <- rvg::dml(ggobj = normal_plot)
    figures_PP <- add_slide(figures_PP,layout = "Title and Content", master = "Office Theme")
    figures_PP <- ph_with(figures_PP, value = slide_title, location = ph_location_type(type = "title"))
    figures_PP <- ph_with(figures_PP, value = vectorized_plot, ph_location_type(type = "body"))
    
    slide_title <- paste0('SingleR Metric: delta next')
    normal_plot <- plotDeltaDistribution(SingleR_predictions, size=NA, show='delta.next') + theme(legend.position = 'none') + plot_annotation(title='SingleR Metric: delta next', theme = theme(plot.title = element_text(hjust = 0.5)))
    vectorized_plot <- rvg::dml(ggobj = normal_plot)
    figures_PP <- add_slide(figures_PP,layout = "Title and Content", master = "Office Theme")
    figures_PP <- ph_with(figures_PP, value = slide_title, location = ph_location_type(type = "title"))
    figures_PP <- ph_with(figures_PP, value = vectorized_plot, ph_location_type(type = "body"))
    
    slide_title <- paste0('Histogram of SingleR Delta Next Values')
    #normal_plot <- hist(SingleR_predictions$delta.next, 100, main='delta.next distribution')
    binwidth_for_histo <- max(SingleR_predictions$delta.next)/100
    histoplot <- ggplot(data.frame(delta.next = SingleR_predictions$delta.next), aes(x=delta.next)) +
      geom_histogram(binwidth=binwidth_for_histo, fill="grey75", color="black", alpha=0.9) +
      theme_classic() + scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
      labs(title = 'Delta.Next Value Distribution')
    #normal_plot <- as.ggplot(normal_plot)
    vectorized_plot <- rvg::dml(ggobj = histoplot)
    figures_PP <- add_slide(figures_PP,layout = "Title and Content", master = "Office Theme")
    figures_PP <- ph_with(figures_PP, value = slide_title, location = ph_location_type(type = "title"))
    figures_PP <- ph_with(figures_PP, value = vectorized_plot, ph_location_type(type = "body"))
    
    save_filename <- paste0(save_dir,'/SingleR_QC_Figs-',Date_Tag,'.pptx')
    print(figures_PP, target = save_filename)
    rm(figures_PP)
    
  }
  
  
  ################################################################
  ################################################################
  ##
  ## Adding SingleR cell type predictions to seurat_data metadata
  
  SingleR_metaprep <- data.frame(labels = SingleR_predictions$labels)
  rownames(SingleR_metaprep) <- rownames(SingleR_predictions)
  
  seurat_data <- AddMetaData(object = seurat_data, metadata = SingleR_metaprep, col.name = 'SingleR_Predictions')
  
  if (save_seurat_SingleR) {
    seurat_output_name <- paste0("seurat_singleR_",data_prefix,"-",Date_Tag)
    saveRDS(seurat_data, paste0(save_dir,"/",seurat_output_name,".rds"))
  }
  
  return(seurat_data)
  
  
  ## Creating SingleR Run Info text file
  SingleR_RunInfo_Filename <- paste0(save_dir,"/",SingleR_output_name,".txt")
  
  sink(file = SingleR_RunInfo_Filename)
  
  print("Run Info for run_SingleR_cell_annotation")
  print("A subfunction in the STEP3_CellTypeIdentification_BCL.R script")
  print("Run on:")
  Sys.time()
  print("-------------------------------------")
  print("Variable and Toggle Selection for:")
  print(paste0(SingleR_output_name,".rds"))
  print("-------------------------------------")
  print("User provided notes on reference dataset used:")
  print(as.character(reference_data_notes))
  print("---------------------------------------------------------------------------------------------")
  print("---------------------------------------------------------------------------------------------")
  print("General Code Option Settings")
  print("----------------------------------------")
  print(paste0("convert_ensembl_to_hgnc = ",convert_ensembl_to_hgnc))
  print(paste0("use_varFeatures_for_SingleR = ",use_varFeatures_for_SingleR))
  print(paste0("n variable features used for SingleR = ",nFeatures_for_SingleR_comp))
  print("----------------------------------------")
  print("Some reference sets contain an unspecified 'native' cell group. This group essentially adds noise to the")
  print("SingleR prediction process. It could theoretically provide additional confidence to cell-type annotations")
  print("made in the presence of this noise, but it can be problematic if the overall strength of SingleR predictions")
  print("is low. The following toggle determines whether this 'native' cell group is filtered out prior to running SingleR")
  print(paste0("remove_native_cells = ",remove_native_cells))
  print("---------------------------------------------------------------------------------------------")
  print("---------------------------------------------------------------------------------------------")
  print("Script Outputs:")
  print(paste0('Output Directory: ',save_dir))
  print("----------------------------------------")
  print(paste0("save_SingleR_predictions_dataframe = ",save_SingleR_predictions_dataframe))
  if (save_SingleR_predictions_dataframe) {
    print(paste0('SingleR dataframe filename: ', SingleR_output_name,'.rds'))
  }
  print("----------------------------------------")
  print(paste0("save_SingleR_QC_figures = ",save_SingleR_QC_figures))
  if (save_SingleR_QC_figures) {
    print(paste0('SingleR QC PowerPoint: ','SingleR_QC_Figs-',Date_Tag,'.pptx'))
  }
  print("----------------------------------------")
  print(paste0("save_seurat_SingleR = ",save_seurat_SingleR))
  if (save_seurat_SingleR) {
    print(paste0('seurat object filename: ',seurat_output_name,'.rds'))
  }
  
  sink(file = NULL)
  
}








##  filter_seurat_using_SingleR_data(seurat_data, SingleR_predictions,
##                                   use_pruned_list=FALSE, delta_med_threshold=-Inf,
##                                   delta_next_threshold=0,
##                                   save_modified_seurat=TRUE, save_dir=getwd())
##
##  This function allows subsetting of the seurat object using metrics from SingleR annotation
##
##  In order to use this function, you will need to have saved the SingleR_predictions data frame
##    created when SingleR was run. You will need to load that file and provide it as one of the
##    inputs to this function.
##
##  seurat_data - The seurat object that you want to subset
##  SingleR_predictions - The SingleR predictions data frame associated with your seurat object.
##      At the moment, you must not have done any previous subsetting of your seurat data. This is
##      something that I can fix in the future.
##      That said, if you used the variable features option when you ran SingleR, then the SingleR
##      QC metrics are likely to change if you re-run your subsetted data. So, if you use the variable
##      features option, then it is probably best to re-run SingleR prior to using this function.
##  use_pruned_list - Logical toggle that determines whether to use SingleR's automated QC prediction
##      thresholding to identify cells that should be subsetted. Note that if the overall prediction
##      scores from SingleR are relatively low, then SingleR's automated QC check will fail. This will
##      result in no cells being identified for removal. In this case, you will be forced to utlize the
##      user defined thresholding options (next two variables)
##  delta_med_threshold - Identify cells for removal using a user-defined delta median threshold value.
##  delta_next_threshold - Identify cell for removal using a user-defined delta next threshold value
##      NOTE: You can threshold using delta_med and/or delta_next. If you only want to use one of the two
##            metrics for thresholding, DO NOT define the other variable (or define it using the default
##            value I've written into the function: -Inf for delta_med or 0 for delta_next)
##  


filter_seurat_using_SingleR_data <- function(seurat_data, SingleR_predictions,
                                             use_pruned_list=FALSE, delta_med_threshold=-Inf,
                                             delta_next_threshold=0, 
                                             save_modified_seurat=TRUE, save_dir=getwd()) {
  
  notes_storage <- seurat_data@meta.data$Notes
  
  modification_check <- grep('Modification',notes_storage)
  if (length(modification_check)>0) {
    modification_num <- length(modification_check)+1
  } else {
    modification_num <- 1
  }
  
  
  if (use_pruned_list) {
    to.remove <- is.na(SingleR_predictions$pruned.labels)
  } else {
    to.remove <- pruneScores(SingleR_predictions, min.diff.med = delta_med_threshold,
                             min.diff.next = delta_next_threshold)
  }
  
  cells_to_remove <- data.frame(cells = rownames(SingleR_predictions), to_remove = to.remove)
  cells_to_remove <- column_to_rownames(cells_to_remove, 'cells')
  seurat_data <- AddMetaData(object = seurat_data, metadata = cells_to_remove, col.name = 'cells_to_remove')
  
  Idents(seurat_data) <- 'cells_to_remove'
  seurat_data <- subset(seurat_data, idents = 'FALSE')
  
  
  ## Updating the Seurat object modification notes
  seurat_data@meta.data$Notes <- notes_storage[1:nrow(seurat_data@meta.data)]
  
  empty_notes_inds <- which(is.na(seurat_data@meta.data$Notes))
  seurat_data@meta.data$Notes[empty_notes_inds[1]] <- '---------------------------'
  seurat_data@meta.data$Notes[empty_notes_inds[2]] <- paste0('Modification ',as.character(modification_num))
  seurat_data@meta.data$Notes[empty_notes_inds[3]] <- 'filter_seurat_using_SingleR_data'
  if (use_pruned_list) {
    seurat_data@meta.data$Notes[empty_notes_inds[4]] <- paste0('use_pruned_list = ',use_pruned_list)
  } else {
    seurat_data@meta.data$Notes[empty_notes_inds[4]] <- paste0('use_pruned_list = ',use_pruned_list)
    seurat_data@meta.data$Notes[empty_notes_inds[5]] <- paste0('delta_med_threshold = ',delta_med_threshold)
    seurat_data@meta.data$Notes[empty_notes_inds[6]] <- paste0('delta_next_threshold = ',delta_next_threshold)
  }
  
  if (save_modified_seurat) {
    new_seurat_name <- paste0(save_dir,'/seurat_data_Modification_',as.character(modification_num),'.rds')
    saveRDS(seurat_data, new_seurat_name)
  }
  
  return(seurat_data)
}

















##  subset_seurat_data(seurat_data, save_modified_seurat=TRUE, save_dir=getwd(),
##                     Ident_ID = 'seurat_clusters', Idents_to_subset, use_invert = TRUE, subset_notes)
##
##  This function will subset your seurat data. The primary purpose of this function is to properly update
##      and maintain the modification history in the Notes column of the seurat metadata.
##
##  seurat_data - The seurat object to be subsetted
##  save_modified_seurat - Toggle to determine whether or not to save the seurat object after subsetting (Default = TRUE)
##  save_dir - The directory where the modified seurat object should be saved (Default = getwd())
##  Ident_ID - The Idents that will be used for subsetting (Default = 'seurat_clusters')
##  Idents_to_subset - Vector of idents that will be provided to Seurat's subset function
##  use_invert - Logical toggle used to set the value for subset's invert option (Default = TRUE)
##  subset_notes - A optional vector of strings containing user notes to be added to the Notes metadata. For example:
##                 subset_notes <- c('Cluster 5 - High mito content', 
##                                   'Cluster 6 - snowman, likely doublets',
##                                   'Cluster 9 - Non-target cell type cluster')


subset_seurat_data <- function(seurat_data, save_modified_seurat=TRUE, save_dir=getwd(),
                               Ident_ID = 'seurat_clusters', Idents_to_subset, use_invert = TRUE, subset_notes) {
  
  
  notes_storage <- seurat_data@meta.data$Notes
  modification_check <- grep('Modification',notes_storage)
  if (length(modification_check)>0) {
    modification_num <- length(modification_check)+1
  } else {
    modification_num <- 1
  }
  
  Idents(seurat_data) <- Ident_ID
  seurat_data <- subset(seurat_data, idents = Idents_to_subset, invert = use_invert)
  
  ## Updating the Seurat object modification notes
  seurat_data@meta.data$Notes <- notes_storage[1:nrow(seurat_data@meta.data)]
  
  empty_notes_inds <- which(is.na(seurat_data@meta.data$Notes))
  seurat_data@meta.data$Notes[empty_notes_inds[1]] <- '---------------------------'
  seurat_data@meta.data$Notes[empty_notes_inds[2]] <- paste0('Modification ',as.character(modification_num))
  seurat_data@meta.data$Notes[empty_notes_inds[3]] <- 'subset_seurat_data'
  seurat_data@meta.data$Notes[empty_notes_inds[4]] <- paste0('Idents = ',Ident_ID)
  seurat_data@meta.data$Notes[empty_notes_inds[5]] <- paste0('use_invert = ',use_invert)
  for (i in 1:length(Idents_to_subset)) {
    seurat_data@meta.data$Notes[empty_notes_inds[5+i]] <- paste0('Subsetted Idents = ',Idents_to_subset[i])
  }
  if (length(subset_notes) > 0) {
    ind_adj <- i+6
    seurat_data@meta.data$Notes[empty_notes_inds[ind_adj]] <- 'User_Defined Subsetting Notes'
    for (i in 1:length(subset_notes)) {
      seurat_data@meta.data$Notes[empty_notes_inds[ind_adj+i]] <- subset_notes[i]
      
    }
  }
  
  if (save_modified_seurat) {
    new_seurat_name <- paste0(save_dir,'/seurat_data_Modification_',as.character(modification_num),'.rds')
    saveRDS(seurat_data, new_seurat_name)
  }
  
  
  return(seurat_data)
}


















##  create_subset_report(seurat_data, pre_sub_figs, post_sub_figs, custom_report_name = 'none', save_dir)
##
##  This function will create a PowerPoint report that includes user-provided figures from before
##      and after seurat object modifications (usually subsetting). You will need to decide which
##      figures you want included in the report and provide them to this function. The latest modification
##      and previous modifications will be read automatically from the Notes metadata of the seurat object
##      and included in the report along with the provided figures.
##
##  seurat_data - Your seurat object with the latest modifications
##  pre_sub_figs - A list of figures from BEFORE the latest modification to your seurat object. This means you
##      will need to have created and saved these figures prior to modifying your seurat data.
##      The idea is that these would be the figures that you used to make the decisions about your latest
##      seurat object modification. If available, this function will use the names(pre_sub_figs) as titles for
##      the slides containing each figure.
##  post_sub_figs - A list of figures from the current seurat object (i.e., the one you provided to this function) 
##      If available, names(pre_sub_figs) will be used as titles for the slides containing each figure.
##  custom_report_name - You can provide a string that will be used as the name for the PowerPoint report. Do not include 
##      '.pptx' to the end of your custom_report_name; the code will add that automatically.
##      If a custom_report_name is not provided, the function will generate a report name in the form of:
##      'STEP3_Mod[current modification number]_Report-[Date tag].pptx'
##  save_dir - The directory where the report should be saved. Do not include "/" at the end of the path to the folder


create_subset_report <- function(seurat_data, pre_sub_figs, post_sub_figs, custom_report_name = 'none', save_dir=getwd()) {
  
  ## Collecting modification information from the STEP3 seurat object metadata
  notes_storage <- seurat_data@meta.data$Notes
  empty_notes_inds <- which(is.na(seurat_data@meta.data$Notes))
  modification_check <- grep('Modification',notes_storage)
  if (length(modification_check)>1) {
    previous_modifications <- modification_check[1:(length(modification_check)-1)]
    current_modification <- modification_check[length(modification_check)]
    
    pre_modification_frame <- data.frame(Mod_Notes = notes_storage[2:(current_modification-1)])
    cur_modification_frame <- data.frame(Mod_Notes = notes_storage[current_modification:(empty_notes_inds[1]-1)])
  } else {
    pre_modification_frame <- data.frame(Mod_Notes = c('Unmodified STEP2 Seurat Object'))
    
    current_modification <- modification_check
    cur_modification_frame <- data.frame(Mod_Notes = notes_storage[current_modification:(empty_notes_inds[1]-1)])
  }
  

  
  #####################
  ##
  ## Creating the PowerPoint object
  
  figures_PP <- read_pptx()
  windowsFonts("Times" = windowsFont("Times"))

  #####################
  ##
  ## Creating the Previous Modification History Slide
  
  slide_title <- paste0('Previous Modification History')
  text_mod <- fp_text_lite(font.size=22, color='black')
  slide_title <- fpar(ftext(slide_title, text_mod), fp_p = fp_par(text.align = 'center'))
  figures_PP <- add_slide(figures_PP,layout = "Two Content", master = "Office Theme")
  #figures_PP <- add_slide(figures_PP,layout = "Title and Content", master = "Office Theme")
  figures_PP <- ph_with(figures_PP, value = slide_title, location = ph_location_type(type = "title"))
  
  ## This is code I may eventually use to make these modification history slides more aesthetically appealing
  #text_mod <- fp_text_lite(font.size=10, color='black')
  #ft_SingleR <- flextable(test)
  #ft_SingleR <- flextable::style(ft_SingleR, pr_t = text_mod)
  #ft_SingleR <- padding(ft_SingleR, padding.top = 0, part = "all")
  #ft_SingleR <- padding(ft_SingleR, padding.bottom = 0, part = "all")
  #ft_SingleR <- autofit(ft_SingleR)
  #figures_PP <- ph_with(figures_PP, value = ft_SingleR, ph_location_type(type='body'))
  
  if (nrow(pre_modification_frame) > 30) {
    slide_content <- paste(pre_modification_frame$Mod_Notes[1:30],collapse ="\n")
    text_mod <- fp_text_lite(font.size=10, color='black')
    slide_content <- fpar(ftext(slide_content, text_mod), fp_p = fp_par(text.align = 'left'))
    figures_PP <- ph_with(figures_PP, value = slide_content, ph_location_left())
    
    slide_content <- paste(pre_modification_frame$Mod_Notes[31:nrow(pre_modification_frame)],collapse ="\n")
    text_mod <- fp_text_lite(font.size=10, color='black')
    slide_content <- fpar(ftext(slide_content, text_mod), fp_p = fp_par(text.align = 'left'))
    figures_PP <- ph_with(figures_PP, value = slide_content, ph_location_right())
    
  } else {
    slide_content <- paste(pre_modification_frame$Mod_Notes,collapse ="\n")
    text_mod <- fp_text_lite(font.size=10, color='black')
    slide_content <- fpar(ftext(slide_content, text_mod), fp_p = fp_par(text.align = 'left'))
    figures_PP <- ph_with(figures_PP, value = slide_content, ph_location_left())
  }
  
  
  #####################
  ##
  ## Adding user provided pre-modification figures
  
  if (length(pre_sub_figs)>0) {
    for (i in 1:length(pre_sub_figs)) {
      
      if (length(names(pre_sub_figs)[i])>0) {
        slide_title <- paste0('Pre-Modification - ',names(pre_sub_figs)[i])
      } else {
        slide_title <- 'Pre-Modification Figure - No Name Provided'
      }
      slide_title <- fpar(ftext(slide_title, text_mod), fp_p = fp_par(text.align = 'center'))
      figures_PP <- add_slide(figures_PP,layout = "Title and Content", master = "Office Theme")
      figures_PP <- ph_with(figures_PP, value = slide_title, location = ph_location_type(type = "title"))
      figures_PP <- ph_with(figures_PP, value = pre_sub_figs[[i]], ph_location_type(type = "body"))
    }
  }

  
  #####################
  ##
  ## Creating the Latest Modification History Slide
  
  slide_title <- paste0('Most Recent Modification')
  text_mod <- fp_text_lite(font.size=22, color='black')
  slide_title <- fpar(ftext(slide_title, text_mod), fp_p = fp_par(text.align = 'center'))
  figures_PP <- add_slide(figures_PP,layout = "Title and Content", master = "Office Theme")
  figures_PP <- ph_with(figures_PP, value = slide_title, location = ph_location_type(type = "title"))
  
  slide_content <- paste(cur_modification_frame$Mod_Notes,collapse ="\n")
  text_mod <- fp_text_lite(font.size=14, color='black')
  slide_content <- fpar(ftext(slide_content, text_mod), fp_p = fp_par(text.align = 'left'))
  figures_PP <- ph_with(figures_PP, value = slide_content, ph_location_type(type = "body"))
  
  
  #####################
  ##
  ## Adding user provided post-modification figures
  
  if (length(post_sub_figs)>0) {
    for (i in 1:length(post_sub_figs)) {
      
      if (length(names(post_sub_figs)[i])>0) {
        slide_title <- paste0('Post-Modification - ',names(post_sub_figs)[i])
      } else {
        slide_title <- 'Post-Modification Figure - No Name Provided'
      }
      slide_title <- fpar(ftext(slide_title, text_mod), fp_p = fp_par(text.align = 'center'))
      figures_PP <- add_slide(figures_PP,layout = "Title and Content", master = "Office Theme")
      figures_PP <- ph_with(figures_PP, value = slide_title, location = ph_location_type(type = "title"))
      figures_PP <- ph_with(figures_PP, value = post_sub_figs[[i]], ph_location_type(type = "body"))
    }
  }
  
  
  
  #########################
  ##
  ##  Saving the PowerPoint report
  
  if (custom_report_name == 'none') {
    thetime <- Sys.time()
    Date_Tag <- paste0(substr(thetime,3,4),substr(thetime,6,7),substr(thetime,9,10),substr(thetime,8,8),
                       substr(thetime,12,13),'h',substr(thetime,15,16),'m')
    
    PP_name <- paste0('STEP3_Mod',as.character(length(modification_check)),'_Report-',Date_Tag,'.pptx')
    save_filename <- paste0(save_dir,'/',PP_name)
  } else {
    save_filename <- paste0(save_dir,'/',custom_report_name)
  }
  
  print(figures_PP, target = save_filename)
  rm(figures_PP)
  
}









##  perform_standard_seurat_workflow(seurat_data, nVariableFeatures = 2000, FindClusterResolution = 0.4,
##                                   variables_to_regress = c('percent_mito','nCount_RNA'),
##                                   use_automated_nPCs_calc = TRUE, PCs_CumulPctVar_Threshold = 80,
##                                   PCs_AbsVar_Threshold = 5, PC_DeltaVar_Threshold = 0.1,
##                                   nPCs_for_cluster_fxns = 15)
##
##  This function is just a convenient way to run a seurat object through the standard Seurat workflow
##      after making modifications. It is set up to use a handful of default values, but you'll probably
##      want to define all of these in your own code and then feed them into this function so they are
##      appropriate for your seurat object.


perform_standard_seurat_workflow <- function(seurat_data, nVariableFeatures = 2000, FindClusterResolution = 0.4,
                                             variables_to_regress = c('percent_mito','nCount_RNA'),
                                             use_automated_nPCs_calc = TRUE, PCs_CumulPctVar_Threshold = 80,
                                             PCs_AbsVar_Threshold = 5, PC_DeltaVar_Threshold = 0.1,
                                             nPCs_for_cluster_fxns = 15) {
  
  
  seurat_data <- seurat_data %>%
    NormalizeData() %>%
    FindVariableFeatures(nfeatures = nVariableFeatures)
  
  ## Removing mitochondrial genes from the list of variable features
  VariableFeatures(seurat_data) <- VariableFeatures(seurat_data)[!grepl("(?i)^mt-", VariableFeatures(seurat_data))]
  
  seurat_data <- seurat_data %>%
    ScaleData(vars.to.regress = variables_to_regress) %>%
    RunPCA()
  
  if (use_automated_nPCs_calc) {
    stdv <- seurat_data[["pca"]]@stdev
    sum.stdv <- sum(seurat_data[["pca"]]@stdev)
    percent.stdv <- (stdv / sum.stdv) * 100
    cumulative <- cumsum(percent.stdv)
    co1 <- which(cumulative > PCs_CumulPctVar_Threshold & percent.stdv < PCs_AbsVar_Threshold)[1]
    co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - 
                         percent.stdv[2:length(percent.stdv)]) > PC_DeltaVar_Threshold), 
                decreasing = T)[1] + 1
    nPCs_for_cluster_fxns <- min(co1, co2)
  }
  
  seurat_data <- seurat_data %>%
    FindNeighbors(reduction = "pca", dims=1:nPCs_for_cluster_fxns) %>%
    FindClusters(resolution = FindClusterResolution) %>%
    RunUMAP(dims = 1:nPCs_for_cluster_fxns, reduction = "pca")
  
  return(seurat_data)
  
}








