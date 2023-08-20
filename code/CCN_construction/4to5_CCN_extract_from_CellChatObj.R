library(dplyr)
library(reticulate)
source_python("../ML/function/DataLoader.py")


## dataset and output path setting
dataset = "Jung"
input_path = file.path("../../data/CCN_construction/4_Commprob_calculated_CellChatObj", dataset)
output_path = file.path("../../data/CCN_construction/5_CCN", dataset)
if(!file.exists(output_path))
    dir.create(output_path, recursive = TRUE)

message(dataset)


## dataset and response info
clinical = nonres_res_partition(dataset)
nonresponder_sample = filter(clinical, response == 0) %>% rownames
responder_sample = filter(clinical, response == 1) %>% rownames

## CellChat object filelist
obj_filelist = list.files(input_path, recursive = FALSE, pattern = ".RData" )
obj_filelist = obj_filelist[ match(rownames(clinical), lapply(strsplit(obj_filelist, "_CellChatObj"), function(x) { x[1] })) ]


## CCN_matrix by pathway
load(file.path(input_path, obj_filelist[1]))
pathway = CellChat@netP$pathways
for(j in 2:length(obj_filelist)){
  load(file.path(input_path, obj_filelist[j]))
  pathway = intersect(pathway, CellChat@netP$pathways)
}
pathway = sort(pathway)


PairedMatrix = list()
for(p in 1:length(pathway)){
    
    for(i in 1:length(obj_filelist)){
        load(file.path(input_path, obj_filelist[i]))
        pat = strsplit(obj_filelist[i], "_CellChatObj")[[1]][1]
        PairedMatrix[[pathway[p]]][[pat]] = list()

        commprob_pathway = CellChat@netP$prob[, , pathway[p]]
        for(m in 1:dim(commprob_pathway)[1]){
            for(n in 1:dim(commprob_pathway)[2]){
                source = rownames(commprob_pathway)[m]
                target = colnames(commprob_pathway)[n]

                PairedMatrix[[pathway[p]]][[pat]][[paste0(source, "_to_", target)]] = commprob_pathway[source, target]
            }
        }
        PairedMatrix[[pathway[p]]][[pat]] = do.call(cbind, PairedMatrix[[pathway[p]]][[pat]]) %>% as.data.frame
    }
    PairedMatrix[[pathway[p]]] = do.call(rbind, PairedMatrix[[pathway[p]]]) %>% as.data.frame
    message(p, " ", pathway[p])

}
CCN_matrix = do.call(cbind, PairedMatrix) %>% as.data.frame
summed_CCN = apply(CCN_matrix, 2, sum)
col = summed_CCN[summed_CCN != 0] %>% names
CCN_matrix = CCN_matrix[, col]


## aggregated and file out
nonresponder_CCN = cbind( sample_id = nonresponder_sample, CCN_matrix[ nonresponder_sample ,] )
responder_CCN = cbind( sample_id = responder_sample, CCN_matrix[ responder_sample ,])
write.table(nonresponder_CCN, file.path(output_path, "nonresponder_CCN.txt"),
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(responder_CCN, file.path(output_path, "responder_CCN.txt"),
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

