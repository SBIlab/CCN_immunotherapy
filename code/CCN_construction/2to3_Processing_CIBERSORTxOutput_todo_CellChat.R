library(dplyr)
library(preprocessCore)


## dataset to process
dataset = "Jung"

## set path
input_path = file.path("../../data/CCN_construction/2_CIBERSORTx_output", dataset)
output_path = file.path("../../data/CCN_construction/3_processed_celltype_gene_expression", dataset)
if(!file.exists(output_path))
    dir.create(output_path)

## get cell type GEP file names
celltypeGEP_filename = list.files(input_path, recursive = FALSE, pattern = "^CIBERSORTxHiRes_NA_")
celltypeGEP_filename = celltypeGEP_filename[ grep(".txt$", celltypeGEP_filename) ]

## load cell type GEP and process
message(paste0("input: ", dataset))
for(i in 1:length(celltypeGEP_filename)){
    exp_df = read.table(file.path(input_path, celltypeGEP_filename[i]), 
                            header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

    # Gene symbol . -> -
    rownames(exp_df) = gsub("\\.", "-", rownames(exp_df))

    # log transformation
    exp_df[is.na(exp_df)] = 0
    exp_df = log(exp_df + 1)

    # quantile normalization
    norm_exp_df = normalize.quantiles(as.matrix(exp_df), copy = FALSE) %>% as.data.frame

    # return processed expression profile
    output_norm_exp_df = cbind( data.frame(GeneSymbol = rownames(norm_exp_df)), norm_exp_df)
    write.table(output_norm_exp_df, file.path(output_path, celltypeGEP_filename[i]),
                sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    message(paste0( celltypeGEP_filename[i], ", processed" ))

}
