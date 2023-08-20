library(dplyr)
library(Matrix)
source("./CellChat-master/R/CellChat_class.R")
source("./CellChat-master/R/modeling_modified.R")
source("./CellChat-master/R/analysis.R")
source("./CellChat-master/R/utilities.R")
source("./CellChat-master/R/database.R")
source("./CellChat-master/R/visualization.R")

## set dataset
dataset = "Jung"

## set path
input_path = file.path("../../data/CCN_construction/3_processed_celltype_gene_expression", dataset)
output_path = file.path("../../data/CCN_construction/4_Commprob_calculated_CellChatObj", dataset)
if(!file.exists(output_path))
  dir.create(output_path, recursive = TRUE)

## load cell type GEP file names and define cell types
deconv_filelist = list.files(path = input_path, recursive = FALSE, pattern = "CIBERSORTxHiRes")
deconv_filelist = deconv_filelist[ grepl(".txt$", deconv_filelist)]
cell_type = lapply(strsplit(deconv_filelist, "_"), function(x) { x[3] }) %>% unlist  

message(paste0("input path: ", input_path))
message(paste0("output path: ", output_path))

## expression data load, define patient name
cell_type_expList = list()
for(i in 1:length(deconv_filelist)){ 
  cell_type_expList[[cell_type[i]]] = read.table(file.path(input_path, deconv_filelist[i]),
                                                  header = TRUE, sep = "\t", check.names = FALSE)
}
patient = colnames(cell_type_expList[[i]])[2:ncol(cell_type_expList[[i]])]
message(paste0("Number of patients: ", length(patient)))

## CellChat DB load
load("/home/ljh621/cell-cell_communication/R_code/CellChat-master/data/CellChatDB.human.rda")
load("/home/ljh621/cell-cell_communication/R_code/CellChat-master/data/PPI.human.rda")
message("CellChat DB load")

## Calculate communication probabilities of each patient and integrate into comm. pathway
for(k in 1:length(patient)){

  ## current patient
  pat = patient[k]
  message("--------------------------------------------------------------------")
  message(paste0(pat, ": ", Sys.time()))

  ## patient cell type GEP extraction  
  pat_exp = data.frame(GeneSymbol = cell_type_expList[[1]]$GeneSymbol)
  for(j in 1:length(cell_type_expList)){
    tmp_df = cell_type_expList[[cell_type[j]]][ , match(c("GeneSymbol", pat), colnames(cell_type_expList[[cell_type[j]]]))]
    colnames(tmp_df)[2] = cell_type[j]
    pat_exp = left_join(pat_exp, tmp_df, by = "GeneSymbol")
  }
  rownames(pat_exp) = pat_exp$GeneSymbol
  pat_exp = pat_exp[, -1] %>% as.matrix
 
  ## cell type as meta data
  meta = data.frame(cell_type = cell_type)
  rownames(meta) = meta$cell_type
  
  ## CellChat object creation
  CellChat = createCellChat(object = pat_exp, meta = meta, group.by = "cell_type")
  
  CellChat@DB = CellChatDB.human
  CellChat = subsetData(CellChat)
  CellChat@var.features[["features"]] = CellChat@data.signaling@Dimnames[[1]]
  CellChat = identifyOverExpressedInteractions(CellChat)
  
  ## communication probability computation and aggregate
  CellChat = computeCommunProb(CellChat)
  if( length(cell_type) <= 4)
    CellChat = computeCommunProbPathway(CellChat, thresh = 1)
  else
    CellChat = computeCommunProbPathway(CellChat)
  CellChat = aggregateNet(CellChat)
  
  ## save CellChat object
  save(CellChat, file = file.path(output_path, paste0(pat, "_CellChatObj.RData")))
  message(paste0("save ", pat, " ", Sys.time()))
  message("--------------------------------------------------------------------")
}



