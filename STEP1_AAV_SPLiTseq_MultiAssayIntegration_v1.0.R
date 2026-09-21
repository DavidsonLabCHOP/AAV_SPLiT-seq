library(ggplot2)
library(reshape2)
library(readr)
library(purrr)
library(tidyverse)
library(data.table)
library(tidyr)
library(dplyr)
library(Seurat)
library(R.utils)
library(tictoc)
library(readxl)


##  The tic() and toc()s throughout the script is code are used to calculate the processing times.
tic('Total Script Run Time')

##  This variable defines the minimum feature counts (genes) threshold used when creating the Seurat
##    object from the mRNA data.
min_feature_threshold <- 100

##  This variable is incorporated into the output file names. 
output_file_tag <- "Retina"

##  This is the directory where output files are saved (the Seurat Object and the AAV counts table)
##      The code will create this directory if it does not already exist
output_dir <- "C:/<Path_to_Your_Directory>/SPLiTseq_Data/seurat"

##  You must provide the full file path to the csv file containing the sample meta data
cell_meta <- "C:/<Path_to_Your_Directory>/SPLiTseq_Data/metadata.csv"

##  This is the base directory that must contain all the necessary data (mRNA tsv.gz counts files and the amplicon counts files)
##      Within this folder you need to create an Excel file that contains the locations of all the data that you would like
##      analyzed. Format requirements for this Excel file are described below.
SPLiTseq_base_dir <- "C:/<Path_to_Your_Directory>/SPLiTseq_Data" 


##  The data_file_id_list variable must contain the full path to the excel file that contains the locations of the data that
##      you want analyzed (I went with full path here to add a little more flexibility)
##  The Excel file must have this basic format: (each set of brackets indicates a different column)
##  [mRNA_File_Path] [Serotype_1_Data_Folder] [Serotype_2_Data_Folder] ... [Serotype_N_Data_Folder]
##  The file itself must contain headers. It doesn't matter what you call the [mRNA_File_Location] column, that gets re-named
##      later in the code. 
##  However, it is important that the column headers for the serotype_data_folders are the names of the corresponding
##      serotype (i.e., if the column contains locations of AAV2 data, then the column header would be AAV2). Those column
##      headers will be used by the script to create labels for the AmpliconAssay and some other outputs.
##  The code will determine the number and identity of serotypes in the dataset based on these column headers. So make sure
##      they are labeled correctly.
##  For the mRNA files, you need to add the path to the tsv.gz file RELATIVE to the SPLiTseq_base_dir. If the files are in
##      the base directory itself, then you just need to add their filenames. If they are located in subfolders within the
##      SPLiTseq_base_dir, then you need to provide the relative path (i.e., if AAV2_mRNA.tsv.gz is located in a subfolder 
##      called AAV2, then the entry would be /AAV2/AAV2_mRNA.tsv.gz or AAV2/AAV2_mRNA.tsv.gz; the code will work with or without
##      a "/" at the start of your filenames and/or paths)
##  For the AAV folders, you need to provide the path to the FOLDER that contains the amplicon counts data files using the
##      same rules described above for the mRNA files. 
##  Lastly, all of this information needs to be in Sheet1 of the workbook.

data_file_id_list <- "C:/<Path_to_Your_Directory>/SPLiTseq_DataFiles.xlsx"
files_list <- data.frame(read_xlsx(data_file_id_list))

mRNA_data_files <- data.frame(mRNA_File = files_list$mRNA_file)

AAV_data_files <- list()
for (i in 2:ncol(files_list)) {
  AAV_data_files[[colnames(files_list)[i]]] <- data.frame(AAV_folder = files_list[,i])
}


##  Creating the output_dir if it does not already exist
if (!dir.exists(output_dir)) {dir.create(output_dir, showWarnings = TRUE, recursive = FALSE, mode = "0777")}


##  Creating a couple of data frames that will be used to track basic stats and processing times.
stat_tracker <- data.frame(matrix(nrow=nrow(mRNA_data_files),ncol=5))
colnames(stat_tracker) <- c('Dataset','min_feature_threshold','mRNA_cells_PreFilter','mRNA_cells_PostFilter','mRNA_cells_Filtered')
for (i in 1:length(AAV_data_files)) {
  AAV_stub <- data.frame(matrix(nrow=nrow(mRNA_data_files),ncol=3))
  colnames(AAV_stub) <- c(paste0(names(AAV_data_files)[i],'_cells_PreFilter'),
                          paste0(names(AAV_data_files)[i],'_cells_PostFilter'),
                          paste0(names(AAV_data_files)[i],'_cells_Filtered'))
  stat_tracker <- cbind(stat_tracker,AAV_stub)
}
                            
time_tracker <- data.frame(matrix(nrow=nrow(mRNA_data_files),ncol=3))
colnames(time_tracker) <- c('Dataset','mRNA_Counts_Loading','Seurat_Object_Creation')
for (i in 1:length(AAV_data_files)) {
  AAV_stub <- data.frame(matrix(nrow=nrow(mRNA_data_files),ncol=2))
  colnames(AAV_stub) <- c(paste0(names(AAV_data_files)[i],'_Table_Creation'),
                          paste0(names(AAV_data_files)[i],'_AmpliconAssay_Creation'))
  time_tracker <- cbind(time_tracker,AAV_stub)
}
time_tracker$Total_Time <- NA


output_files_tracker <- c()

for (i in 1:nrow(mRNA_data_files)) {
  
  tic('Total File Processing Time')
  
  ## This code attempts to pull the sample ID from the mRNA tsv.gz counts file name.
  ##    NOTE: this is the portion of the code that is most likely to need editing in the future.
  ##    Currently, it is looking for a sample ID in the form of [One Letter][One Number]. 
  ##    If we change the format of the sample ID then we will need to update the grep() 
  ##    call accordingly. Failure to do so will either cause the code to crash or give you
  ##    an NA in your filenames (which will end up causing subsequent sublibraries from the same
  ##    serotype to over-write the previous sublibrary's data)
  name_test <- str_split(mRNA_data_files$mRNA_File[i],pattern=c('/'))
  name_test <- str_split(name_test[[1]][length(name_test[[1]])],pattern=c('_'))
  data_set_id <- name_test[[1]][grep('^[A-Z][0-9]',name_test[[1]])]
  

  ##  The output file name is constructed here. It currently has the format:
  ##      [output_file_tag]_[data_set_id]
  output_filename <- paste0(output_file_tag,'_',data_set_id)
  stat_tracker$Dataset[i] <- output_filename
  time_tracker$Dataset[i] <- output_filename
  
  stat_tracker$min_feature_threshold[i] <- min_feature_threshold
  
  print(paste0("Processing ",mRNA_data_files$mRNA_File[i],". Save file name is ", paste0(output_filename, "_seurat.rds")))
  

  
  if (substr(mRNA_data_files$mRNA_File[i],1,1) == '/') {
    current_mRNA_file <- paste0(SPLiTseq_base_dir,mRNA_data_files$mRNA_File[i])
  } else {
    current_mRNA_file <- paste0(SPLiTseq_base_dir,'/',mRNA_data_files$mRNA_File[i])
  }


  Sys.setenv("VROOM_CONNECTION_SIZE" = 2000000000) 

  
  ## Load in mRNA counts file data
  ##    I changed the file import method. The previous code used read_delim. I've changed it to use
  ##    fread instead because it is over 10x faster than read_delim
  tic('mRNA counts loading')
  counts <- as_tibble(data.table::fread(current_mRNA_file))
  counts <- column_to_rownames(counts, var="gene")
  time_check <- toc()
  time_tracker$mRNA_Counts_Loading[i] <- time_check$callback_msg

  print('mRNA data loaded')

  
  ## Filling in the stats tracker. I know it is gratuitous to create the mRNA_cell_count_prefilter
  ##    variable, but I did it this way so it is easier to see where the info is coming from.
  mRNA_cell_count_prefilter <- ncol(counts)
  stat_tracker$mRNA_cells_PreFilter[i] <- mRNA_cell_count_prefilter
  
  # Read in metadata
  metadata <- read_csv(cell_meta, show_col_types = FALSE)
  
  # Create Seurat Object with FeatureCounts threshold implemented
  tic('Seurat Object creation')
  seurat <- CreateSeuratObject(counts, min.features=min_feature_threshold) ## add in min cells
  rm(counts)
  time_check <- toc()
  time_tracker$Seurat_Object_Creation[i] <- time_check$callback_msg
  
  print('Seurat object created')
  
  ##  This meta data updating code is fast so I didn't bother tracking it's processing time.
  CellsMeta <- seurat@meta.data
  
  # Create a vector the first 8 nucleotides of each cell barocde
  Barcode <- substr(rownames(CellsMeta), 1, 8)
  
  #table(Barcode == substr(rownames(CellsMeta), 1, 8))
  
  # Add cell round 1 barcode to seurat object metadata
  seurat<- AddMetaData(object = seurat,metadata = Barcode, col.name = "Barcode")
  
  #table(seurat$Barcode == substr(rownames(CellsMeta), 1, 8))
  #table(substr(rownames(CellsMeta), 1, 8) == seurat$Barcode)
  
  # Merge cell metadata via cell round 1 barcode
  seurat@meta.data <- merge(seurat@meta.data, metadata, by = "Barcode", all = T)
  seurat@meta.data <- seurat@meta.data[complete.cases(seurat@meta.data [,c("orig.ident")]),]
  rownames(seurat@meta.data)<-rownames(CellsMeta)
  
  seurat$orig.ident <- output_filename
  
  
  ## Pulling out the mRNA cellIDs that remain after applying the FeatureCounts threshold
  mRNA_cells <- Cells(seurat)
  
  ## Filling in the stats tracker. Again, it's gratuitous to create the mRNA_cell_count_postfilter
  ##    variable, but I did it this way so it is easier to see where the info is coming from.
  mRNA_cell_count_postfilter <- length(mRNA_cells)
  stat_tracker$mRNA_cells_PostFilter[i] <- mRNA_cell_count_postfilter
  stat_tracker$mRNA_cells_Filtered[i] <- mRNA_cell_count_prefilter - mRNA_cell_count_postfilter
  
  
  ##  Loading in the Amplicon Data
  ##  The code will create and add an AmpliconAssay for all serotypes associated with a given
  ##      mRNA data file.
  
  for (x in 1:length(AAV_data_files)) {
    
    serotype <- names(AAV_data_files)[x]
    
    ## These if/else statements check whether or not the entries in the .txt files begin with "/"
    if (substr(AAV_data_files[[x]][i,1],1,1) == '/') {
      current_aav_dir <- paste0(SPLiTseq_base_dir,AAV_data_files[[x]][i,1])
    } else {
      current_aav_dir <- paste0(SPLiTseq_base_dir,'/',AAV_data_files[[x]][i,1])
    }
    
    ## Creation of list of amplicon cellIDs from the contents of the amplicon folder.
    files <-list.files(path = current_aav_dir) %>%
      as_tibble() %>%
      dplyr::rename(name = value)
    
    ## Pulling pre filter AAV cell count for the stat_tracker from the length of the AAV files list
    aav_cell_count_prefilter <- nrow(files)
    
    ## Speeding up the code by pre-filtering the list of amplicon data files to remove
    ##    any cellIDs that have already been filtered out by the FeatureCounts threshold
    ##    applied when creating the seurat object
    files <- as_tibble(files$name[files$name %in% mRNA_cells])
    colnames(files) <- c('name')
    
    ## Filling in the stats tracker. Again, gratuitous variables used for ease of identifying source of info.
    ##    The files variable has been reduced so that it only contains AAV cell ids found in the filtered
    ##    mRNA cells list, so it now provides the post filter AAV cell count.
    aav_cell_count_postfilter <- nrow(files)
    stat_tracker[i,6+(x-1)*3] <- aav_cell_count_prefilter
    stat_tracker[i,7+(x-1)*3] <- aav_cell_count_postfilter
    stat_tracker[i,8+(x-1)*3] <- aav_cell_count_prefilter - aav_cell_count_postfilter
    
    print(paste0('Begining ',serotype,' counts file loading'))
    
    
    tic(paste0(serotype,' data loading and results table creation'))
    results_table <- as_tibble(fread(paste(current_aav_dir,files$name[1],sep="/")))
    for (y in 2:nrow(files)) {
      results_table <- full_join(results_table,as_tibble(fread(paste(current_aav_dir,files$name[y],sep="/"))), by="AAsequence")
      results_table <- results_table[!duplicated(results_table$AAsequence),]
    }
    
    ## This is.na() step actually takes a while for these large tables. It is possible to speed all of this up a lot if we
    ##    were to change the output of the python scripts so that the amplicon counts files contained entries for
    ##    all the oPool barcodes and just had 0's for barcodes that were not detected.
    results_table[is.na(results_table)] <- 0
    time_check <- toc()
    time_tracker[i,4+(x-1)*2] <- time_check$callback_msg
    
    
    print(paste0(serotype,' Counts File loading complete, begining processing to create AmpliconAssay'))
    
    tic(paste0('Creating data frame for ',serotype,' Seurat Amplicon Assay'))
    
    ## The following code provides a massive speed increase over the previous method. It takes seconds
    ##    to perform the function that previously took so long it would crash my computer.
    non_aav_cells <- data.frame(matrix(nrow=nrow(results_table),ncol=length(mRNA_cells)))
    colnames(non_aav_cells) <- mRNA_cells
    non_aav_cells[,] <- 0
    
    dup_cellids <- colnames(non_aav_cells) %in% colnames(results_table)
    non_aav_cells <- non_aav_cells[,!dup_cellids]
    
    results_table_allcells <- cbind(results_table,non_aav_cells)
    rm(results_table)
    results_table_allcells <- column_to_rownames(results_table_allcells, var = 'AAsequence')
    results_table_allcells <- results_table_allcells[,colnames(results_table_allcells) %in% mRNA_cells]
    
    ##  The code I've enclosed in the following if statements take significant time to execute on these huge
    ##      data frames and it is usually unnecessary.
    ##  If these lines of code definitely needed to be run, then it would be faster to not run the checks first.
    ##      However, because they are usually unnecessary, it ends up saving time to test the data first using
    ##      these much faster checks.
    if (sum(is.na(results_table_allcells)=='TRUE') > 0) {
      results_table_allcells[is.na(results_table_allcells)] <- 0
    }
    
    if (sum((results_table_allcells-floor(results_table_allcells)==0)=='FALSE') > 0) {
      results_table_allcells <- mutate_all(results_table_allcells, as.integer)
    }
    
    time_check <- toc()
    time_tracker[i,5+(x-1)*2] <- time_check$callback_msg
    
    current_dir <- getwd()
    setwd(output_dir)
    saveRDS(results_table_allcells, file = paste0(output_filename,"_",serotype,"_counts_allcells.rds"))
    output_files_tracker <- c(output_files_tracker, paste0(output_filename,"_",serotype,"_counts_allcells.rds"))
    setwd(current_dir)
    
    amplicon_assay <- CreateAssayObject(counts = results_table_allcells, check.matrix = TRUE)
    rm(results_table_allcells)
    
    assay_name <- paste0(serotype,'_AMPLICON')
    seurat[[assay_name]] <- amplicon_assay
    
    rm(amplicon_assay)
    
  }
  
  #seurat@assays
  #head(colnames(seurat))
  #table(seurat$Region)
  #head(seurat@meta.data)
  
  setwd(output_dir)
  saveRDS(seurat, file = paste0(output_filename, "_seurat.rds"))
  output_files_tracker <- c(output_files_tracker, paste0(output_filename, "_seurat.rds"))
  setwd(current_dir)
  
  rm(seurat, CellsMeta, metadata)
  
  time_check <- toc()
  time_tracker$Total_Time[i] <- time_check$callback_msg
  
}

loop_runtime <- toc()


####################################################################
##
##  Creation and Saving of the Output Stats Report
##
####################################################################

thetime <- Sys.time()
Date_Tag <- paste0(substr(thetime,3,4),substr(thetime,6,7),substr(thetime,9,10),substr(thetime,8,8),
                   substr(thetime,12,13),'h',substr(thetime,15,16),'m-',"EST")

Stat_Output_Filename <- paste0("RunInfoAndStats-STEP1-",output_file_tag,"-",Date_Tag,".txt")
setwd(output_dir)

sink(file = Stat_Output_Filename)

print("STEP1_AAV_SPLiTseq_MultiAssayIntegration_v1.R")
print("Run on:")
Sys.time()
print("-------------------------------------")
print(paste("SPLiTseq Data Root Directory:",SPLiTseq_base_dir))
print(paste("Cell MetaData File:",cell_meta))
print("---------------------------------------------------------------------------------------------")
print("Stats of Processed Files:")
print("---------------------")
for (i in 1:nrow(stat_tracker)) {
  print(paste("Dataset Name:",as.character(stat_tracker$Dataset[i])))
  print(paste("mRNA Data File: ",as.character(mRNA_data_files$mRNA_File[i])))
  print(paste("Minimum Gene (Feature) Count Threshold:",as.character(stat_tracker$min_feature_threshold[i])))
  print(paste("Pre-Filter mRNA Cell Count:",as.character(stat_tracker$mRNA_cells_PreFilter[i])))
  print(paste("Post-Filter mRNA Cell Count:",as.character(stat_tracker$mRNA_cells_PostFilter[i])))
  print(paste("Number of mRNA Cells Filtered:",as.character(stat_tracker$mRNA_cells_Filtered[i])))
  for (x in 1:length(AAV_data_files)) {
    print("---------------------")
    print(paste("Amplicon Data Folder:",as.character(AAV_data_files[[x]][i,1])))
    print(paste0("Pre-Filter ",names(AAV_data_files)[x]," Amplicon Cell Count: ",as.character(stat_tracker[i,6+(x-1)*3])))
    print(paste0("Post-Filter ",names(AAV_data_files)[x]," Amplicon Cell Count: ",as.character(stat_tracker[i,7+(x-1)*3])))
    print(paste0("Number of ",names(AAV_data_files)[x]," Amplicon Cells Filtered: ",as.character(stat_tracker[i,8+(x-1)*3])))
  }
  print("---------------------------------------------------------------------------------------------")
  print(paste("Processing Time Stats for",as.character(time_tracker$Dataset[i])))
  print(as.character(time_tracker$mRNA_Counts_Loading[i]))
  print(as.character(time_tracker$Seurat_Object_Creation[i]))
  for (x in 1:length(AAV_data_files)) {
    print(as.character(time_tracker[i,4+(x-1)*2]))
    print(as.character(time_tracker[i,5+(x-1)*2]))
  }
  print(as.character(time_tracker$Total_Time[i]))
  print("---------------------------------------------------------------------------------------------")
  print("---------------------------------------------------------------------------------------------")
}

print(loop_runtime$callback_msg)
print("---------------------------------------------------------------------------------------------")
print("---------------------------------------------------------------------------------------------")
print("Data Files Output:")
print("-----------------------------")
print(paste0("Output Directory: ",output_dir))
for (i in 1:length(output_files_tracker)) {
  print(as.character(output_files_tracker[i]))
}

sink(file = NULL)


## Add input file names (and directory) to the RunInfoAndStats output file.
## For the above, may be easier to read if the root directory is output and then just the file
##    names off that root directory
## Add output file names (and directory)
## Add option to choose whether to save the amplicon counts matrices.



