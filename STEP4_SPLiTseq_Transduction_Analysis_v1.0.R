## This code contains functions designed to analyze the AAV content of AAV_SPLiTseq experimental data

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
  library(DoubletFinder)
  library(presto)
  
  ## Other
  library(officer)
  library(rvg)
  library(writexl)
  library(ggplotify)
  
})




output_dir <- "C:/<Path_To_Directory_Where_Outputs_Will_Be_Saved>"


library(writexl)
library(ggplotify)


seurat_data <- LoadSeuratRds("<Path_to_Your_Seurat_Object_Processed_with_STEP1_and_STEP2_Scripts>")
reduction_to_plot <- 'harmony.umap'
#reduction_to_plot <- 'pca'
#reduction_to_plot <- 'harmony'

table(seurat_data@meta.data$SingleR_predictions)
Seurat::Assays(seurat_data)


#############################
##
##  Code to take a general look at distribution of AAV transduced cells vs non-transduced cells
##

serotype <- '<String that matches the name of your serotype label>'
transduced_cells <- which(seurat_data@meta.data$AAV2_transduction == TRUE)
untransduced_cells <- which(seurat_data@meta.data$AAV2_transduction == FALSE)  

cells_to_highlight <- list()
cells_to_highlight[['Transduced']] <- transduced_cells
#cells_to_highlight[['Not_Transduced']] <- untransduced_cells
#levels(cells_to_highlight) <- c('Transduced','Not_Transduced')

Idents(seurat_data) <- "cell_type"
#transduction_plot <- DimPlot(object = seurat_data, reduction = reduction_to_plot, 
transduction_plot <- DimPlot(object = seurat_data, reduction = reduction_to_plot, 
                             group.by = c("cell_type"),
                             #cols = rep('grey50',length(levels(Idents(seurat_data)))),
                             cells.highlight = cells_to_highlight,
                             #cols.highlight = c('red','grey'),
                             sizes.highlight = c(0.5),
                             label = TRUE) +
plot_annotation(title = paste0(serotype," Transduction Status")) & theme(plot.title = element_text(hjust = 0.5, size=18))

transduction_plot







#############################
##
##  Code to make a DimPlot with individual Barcodes highlighted
##

Assay_highlight_lists <- list()
AAV_labels_lists <- list()
AAV_BCs_list <- list()

assay_names <- c()
for (x in 2:length(Seurat::Assays(seurat_data))) {
  assay_name <- Seurat::Assays(object = seurat_data)[x]
  assay_names <- c(assay_names, assay_name)
  AAV_BCs_list[[assay_name]] <- rownames(seurat_data[[assay_name]])
  
}


for (i in 1:length(assay_names)) {
  
  AAV_labels <- c()
  AAV_highlight_list <- list()
  
  for (x in 1:length(AAV_BCs_list[[assay_names[i]]])) {
    
    DefaultAssay(seurat_data) <- assay_names[i]
    indiv_bc_data <- FetchData(object = seurat_data, vars = c(AAV_BCs_list[[assay_names[i]]][x]), layer = "counts")
    bc_meta_label <- paste0(str_split(assay_names[i],pattern ="_")[[1]][1],'_',colnames(indiv_bc_data))
    indiv_bc_data$temp <- indiv_bc_data[,1] > 0
    indiv_bc_data[,1] <- NULL
    colnames(indiv_bc_data) <- bc_meta_label
    #seurat_data <- AddMetaData(seurat_data, indiv_bc_data)
    AAV_labels <- c(AAV_labels, bc_meta_label)
    AAV_highlight_list[[AAV_BCs_list[[assay_names[i]]][x]]] <- which(indiv_bc_data == TRUE)
    
  }
  
  Assay_highlight_lists[[assay_names[i]]] <- AAV_highlight_list
  AAV_labels_lists[[assay_names[i]]] <- AAV_labels
  
}


DefaultAssay(seurat_data) <- 'RNA'

index_set_lists <- list()
group_size <- 5

for (i in 1:length(assay_names)) {
  index_check <- 1
  index_list <- list()
  while (index_check < length(AAV_labels_lists[[assay_names[i]]])) {
    index_set <- c(index_check, index_check+group_size-1)
    index_check <- index_check + group_size
    index_list[[length(index_list)+1]] <- index_set
  }
  if (index_check >= length(AAV_labels_lists[[assay_names[i]]])) {
    if (index_check == length(AAV_labels_lists[[assay_names[i]]])) {
      index_set <- c(index_check,index_check)
      index_list[[length(index_list)+1]] <- index_set
    } else {
      index_set <- c(index_check,length(AAV_labels_lists[[assay_names[i]]]))
      index_list[[length(index_list)+1]] <- index_set
    }
  }
  index_set_lists[[assay_names[i]]] <- index_list
}



highlight_colors <- c('#3CB22D',
                      '#FED976',
                      '#CC5252',
                      '#4EB3D3',
                      '#A852CC')

transduction_plot_lists <- list()

for (i in 1:length(assay_names)) {
  
  transduction_plots <- list()
  serotype <- str_split(assay_names[i],pattern ="_")[[1]][1]
  index_sets <- index_set_lists[[assay_names[i]]]
  AAV_highlight_list <- Assay_highlight_lists[[assay_names[i]]]
  AAV_labels <- AAV_labels_lists[[assay_names[i]]]
  
  metadata_index <- paste0(serotype,'_transduction')
  cells_to_highlight <- which(seurat_data@meta.data[,metadata_index] == TRUE)
  
  Idents(seurat_data) <- "cell_type"
  transduction_plot <- DimPlot_scCustom(seurat_object = seurat_data, reduction = reduction_to_plot, 
                                        group.by = c("cell_type"),
                                        #colors_use = rep('grey90',length(levels(Idents(seurat_data)))),
                                        cells.highlight = list(cells_to_highlight),
                                        #cols.highlight = highlight_colors,
                                        sizes.highlight = c(0.5),
                                        label = TRUE,
                                        repel = TRUE) +
    labs(title = element_blank()) +
    scale_color_manual(labels = c('Untransduced','Transduced'), values = c('grey80',highlight_colors[3])) +
    plot_annotation(title = paste0(serotype," Transduction Status")) & theme(plot.title = element_text(hjust = 0.5, size=18))
  transduction_plots[[1]] <- transduction_plot
  
  for (x in 1:length(index_set_lists[[assay_names[i]]])) {
    
    index_set <- index_sets[[x]]
    cells_to_highlight <- AAV_highlight_list[index_set[1]:index_set[2]]
    highlight_color_set <- highlight_colors[1:(index_set[2]-index_set[1]+1)]
    
    Idents(seurat_data) <- "cell_type"
    transduction_plot <- DimPlot_scCustom(seurat_object = seurat_data, reduction = reduction_to_plot, 
                                          group.by = c("cell_type"),
                                          #colors_use = rep('grey90',length(levels(Idents(seurat_data)))),
                                          cells.highlight = cells_to_highlight,
                                          #cols.highlight = highlight_colors,
                                          #sizes.highlight = c(0.5),
                                          label = TRUE,
                                          repel = TRUE) +
      labs(title = element_blank()) +
      scale_color_manual(labels = c('Untransduced',names(AAV_highlight_list[index_set[1]:index_set[2]])), values = c('grey80',highlight_color_set)) +
      plot_annotation(title = paste0(serotype," Transduction Status")) & theme(plot.title = element_text(hjust = 0.5, size=18))
    
    transduction_plots[[x+1]] <- transduction_plot
    
  }

  transduction_plot_lists[[assay_names[i]]] <- transduction_plots
  
}


if (save_figures) {
  
  figures_PP <- read_pptx()
  #windowsFonts("Times" = windowsFont("Times"))
  #text_mod <- fp_text_lite(font.size=22, color='black')
  
  for (i in 1:length(transduction_plot_lists)) {
    
    transduction_plots <- transduction_plot_lists[[i]]
    
    for (x in 1:length(transduction_plots)) {
      figures_PP <- add_slide(figures_PP,layout = "Title and Content", master = "Office Theme")
      figures_PP <- ph_with(figures_PP, value = transduction_plots[[x]], ph_location_fullsize())
    }
  }
  
  save_filename <- paste0(output_dir,'/AAV_Transduction_UMAPs.pptx')
  print(figures_PP, target = save_filename)
  rm(figures_PP)
  
}




########################################################
##
##  This code will automatically calculate a number of basic AAV transduction
##      stats. 
##

#seurat_data@meta.data$cell_type <- seurat_data@meta.data$SingleR_predictions


metadata <- seurat_data@meta.data

cell_type_count_threshold <- 100

cell_type_count_table <- as.data.frame(table(seurat_data@meta.data$cell_type))
cell_type_count_table <- cell_type_count_table[cell_type_count_table$Freq > cell_type_count_threshold,]

cell_type_list <- cell_type_count_table$Var1


AAV_BCs_list <- list()
assay_names <- c()
for (x in 2:length(Seurat::Assays(seurat_data))) {
  assay_name <- Seurat::Assays(object = seurat_data)[x]
  assay_names <- c(assay_names, assay_name)
  AAV_BCs_list[[assay_name]] <- rownames(seurat_data[[assay_name]])

}



sum_UMIs_stats <- list()
nTransduced_stats <- list()
pctTransduced_stats <- list()
avg_transd_UMIs_stats <- list()
median_transd_UMIs_stats <- list()


for (i in 1:length(assay_names)) {
  
  aav_data_summary <- data.frame(matrix(nrow=length(AAV_BCs_list[[i]]), ncol=length(cell_type_list)))
  rownames(aav_data_summary) <- AAV_BCs_list[[i]]
  colnames(aav_data_summary) <- cell_type_list
  
  sum_UMIs_stats[[assay_names[i]]] <- aav_data_summary
  nTransduced_stats[[assay_names[i]]] <- aav_data_summary
  pctTransduced_stats[[assay_names[i]]] <- aav_data_summary
  avg_transd_UMIs_stats[[assay_names[i]]] <- aav_data_summary
  median_transd_UMIs_stats[[assay_names[i]]] <- aav_data_summary
  
}


Idents(seurat_data) <- 'cell_type'

for (i in 1:length(cell_type_list)) {
  seurat_subset <- subset(seurat_data, idents = cell_type_list[i])
  
  for (x in 1:(length(Seurat::Assays(seurat_data))-1)) {
    assay_name <- Seurat::Assays(object = seurat_data)[x+1]
    
    for (y in 1:length(AAV_BCs_list[[assay_name]])) {
      
      sum_UMIs_stats[[assay_name]][y,i] <- sum(seurat_subset[[assay_name]]$counts[y,])
      
      nTransduced_stats[[assay_name]][y,i] <- sum(seurat_subset[[assay_name]]$counts[y,]>0)
      pctTransduced_stats[[assay_name]][y,i] <- sum(seurat_subset[[assay_name]]$counts[y,]>0) / length(Cells(seurat_subset))

      transduced_umis <- seurat_subset[[assay_name]]$counts[y,seurat_subset[[assay_name]]$counts[y,]>0]
      if (length(transduced_umis)>0) {
        avg_transd_UMIs_stats[[assay_name]][y,i] <- mean(transduced_umis)
        median_transd_UMIs_stats[[assay_name]][y,i] <- median(transduced_umis)
      } else {
        avg_transd_UMIs_stats[[assay_name]][y,i] <- 0
        median_transd_UMIs_stats[[assay_name]][y,i] <- 0
      }
      
    }
  }
}



format_excel_output <- function(data) {
  data <- rownames_to_column(data, var = 'Peptide')
}


sum_UMIs_stats <- lapply(sum_UMIs_stats, format_excel_output)
nTransduced_stats <- lapply(nTransduced_stats, format_excel_output)
pctTransduced_stats <- lapply(pctTransduced_stats, format_excel_output)
avg_transd_UMIs_stats <- lapply(avg_transd_UMIs_stats, format_excel_output)
median_transd_UMIs_stats <- lapply(median_transd_UMIs_stats, format_excel_output)








####################################################################
##
##  Save basic stats to Excel files
##

excel_output_save_dir <- "<Path_to_Your_Directory>"
Save_XLS <- paste0(excel_output_save_dir,"/sum_UMIs_stats.xlsx")
write_xlsx(sum_UMIs_stats, Save_XLS)

Save_XLS <- paste0(excel_output_save_dir,"/nTransduced_stats.xlsx")
write_xlsx(nTransduced_stats, Save_XLS)

Save_XLS <- paste0(excel_output_save_dir,"/pctTransduced_stats.xlsx")
write_xlsx(pctTransduced_stats, Save_XLS)

Save_XLS <- paste0(excel_output_save_dir,"/avg_transd_UMIs_stats.xlsx")
write_xlsx(avg_transd_UMIs_stats, Save_XLS)

Save_XLS <- paste0(excel_output_save_dir,"/median_transd_UMIs_stats.xlsx")
write_xlsx(median_transd_UMIs_stats, Save_XLS)






####################################################################
##
##  Code to create heatmaps from the basic stats data
##


c.pal <- colorRampPalette(c("white","lavenderblush","thistle1","darkorchid","maroon4"))(10000)

plot_heatmap_justhits <- function(Heat_input) {
  if (nrow(Heat_input>20)) {
    rowfontsize = 8
  } else {
    rowfontsize = 10
  }
  
  Heat_input[Heat_input==0] <- NA
  
  #Heat_input2 <- Heat_input[,4:ncol(Heat_input)]
  
  #find_data <- grep('AVG',colnames(Heat_input))
  
  Heat_input$AVG <- rowMeans(Heat_input[2:ncol(Heat_input)])
  Heat_input <- dplyr::arrange(Heat_input, desc(AVG))
  
  rownames(Heat_input) <- Heat_input$Peptide
  output_graph <- pheatmap(Heat_input[,2:(ncol(Heat_input)-1)],        
                           labels_row = Heat_input$Peptide,
                           cluster_cols = FALSE,
                           color = c.pal,
                           cluster_rows = FALSE,
                           dendrogram="none",
                           trace="none",
                           scale = "none",
                           fontsize = rowfontsize,
                           margins = c(2, 2),
                           border_color=NA,
                           na_col = "white",
                           angle_col = 315,
                           #gaps_row = 27,
                           fontsize_row = rowfontsize)
  return(output_graph)  
}



library(pheatmap)
#test <- plot_heatmap_justhits(pctTransduced_stats[[1]])
#test <- plot_heatmap_justhits(pctTransduced_stats[[2]])

#test <- plot_heatmap_justhits(avg_transd_UMIs_stats[[1]])
#test <- plot_heatmap_justhits(avg_transd_UMIs_stats[[2]])

#test <- plot_heatmap_justhits(median_transd_UMIs_stats[[1]])
#test <- plot_heatmap_justhits(median_transd_UMIs_stats[[2]])




c.pal <- colorRampPalette(c("khaki","yellow",'darkgoldenrod1',"orange",'red','red2',"red3","firebrick","maroon4"))(10000)

plot_heatmap_alt <- function(Heat_input, scale_method, serotype = 'AAV', metric = '') {
  if (nrow(Heat_input>20)) {
    rowfontsize = 8
  } else {
    rowfontsize = 10
  }
  
  Heat_input[Heat_input==0] <- NA
  
  heatmap_title <- paste0(serotype,' - ',metric,'\nnormalization scaling = ',scale_method)
  
  #Heat_input2 <- Heat_input[,4:ncol(Heat_input)]
  
  #find_data <- grep('AVG',colnames(Heat_input))
  
  Heat_input$AVG <- rowMeans(Heat_input[2:ncol(Heat_input)])
  Heat_input <- dplyr::arrange(Heat_input, desc(AVG))
  
  rownames(Heat_input) <- Heat_input$Peptide
  output_graph <- pheatmap(Heat_input[,2:(ncol(Heat_input)-1)],        
                           labels_row = Heat_input$Peptide,
                           cluster_cols = FALSE,
                           #color = c.pal,
                           cluster_rows = FALSE,
                           dendrogram="none",
                           trace="none",
                           main = heatmap_title,
                           scale = scale_method,
                           fontsize = rowfontsize,
                           margins = c(2, 2),
                           border_color=NA,
                           #na_col = "darkblue",
                           na_col = 'navyblue',
                           angle_col = 315,
                           #gaps_row = 27,
                           fontsize_row = rowfontsize)
  
  output_graph <- as.ggplot(output_graph)
  
  return(output_graph)  
}


#nTrans_plot <- plot_heatmap_alt(pctTransduced_stats[[i]], 'none', serotype, 'n Transduced Cells')
#nTrans_plot



save_transduction_heatmaps <- function(nTransduced_stats, pctTransduced_stats) {
  
  figures_PP <- read_pptx()
  #windowsFonts("Times" = windowsFont("Times"))
  #text_mod <- fp_text_lite(font.size=22, color='black')
  
  for (i in 1:length(nTransduced_stats)) {
    
    serotype <- str_split(names(nTransduced_stats)[i], pattern = '_')[[1]][1]
    
    nTrans_plot <- plot_heatmap_alt(nTransduced_stats[[i]], 'column', serotype, 'n Transduced Cells')
    vectorized_plot <- rvg::dml(ggobj = nTrans_plot)
    figures_PP <- add_slide(figures_PP,layout = "Title and Content", master = "Office Theme")
    figures_PP <- ph_with(figures_PP, value = vectorized_plot, ph_location_fullsize())
    
    pctTrans_plot <- plot_heatmap_alt(pctTransduced_stats[[i]], 'none', serotype, '% Transduced Cells')
    vectorized_plot <- rvg::dml(ggobj = pctTrans_plot)
    figures_PP <- add_slide(figures_PP,layout = "Title and Content", master = "Office Theme")
    figures_PP <- ph_with(figures_PP, value = vectorized_plot, ph_location_fullsize())
    
  }
  
  save_filename <- paste0(output_dir,'/AAV_Transduction_Heatmaps.pptx')
  print(figures_PP, target = save_filename)
  rm(figures_PP)
}



if (save_figures) {
  save_transduction_heatmaps(nTransduced_stats, pctTransduced_stats)
}



#test <- plot_heatmap_alt(nTransduced_stats[[1]], 'column')
#test <- plot_heatmap_alt(nTransduced_stats[[2]], 'column')

#test <- plot_heatmap_alt(pctTransduced_stats[[1]], 'none')
#test <- plot_heatmap_alt(pctTransduced_stats[[2]], 'none')

#test <- plot_heatmap_alt(avg_transd_UMIs_stats[[1]], 'column')
#test <- plot_heatmap_alt(avg_transd_UMIs_stats[[2]], 'column')

#test <- plot_heatmap_alt(median_transd_UMIs_stats[[1]], 'column')
#test <- plot_heatmap_alt(median_transd_UMIs_stats[[2]], 'column')

#test <- plot_heatmap_alt(sum_UMIs_stats[[1]], 'column')
#test <- plot_heatmap_alt(sum_UMIs_stats[[2]], 'column')







