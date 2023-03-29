import pandas as pd
import numpy as np
import glob
import os
import itertools
exec(open("./function/DataLoader.py").read())
exec(open("./function/Control_DataLoader.py").read())
exec(open("./function/ML.py").read())

## options
dataset = "VanAllen"  # VanAllen Liu Gide Hugo Jung Kim Mariathasan
result_dir = "../../result/ML/" + dataset + "/"
if not os.path.exists(result_dir):
    os.mkdir(result_dir)
    
gene_num = 3000
qval = 0.05


## load data
perf_list, result_list = [], []
# geneset_name, geneset = load_controlGeneset()

CCN = load_CCN(dataset, gene_num=gene_num, p_cut_off=qval, standardize="StandardScaler", NetworkSelection=True)
clinical = nonres_res_partition(dataset).loc[list(CCN.index), :]
# bulkseq_data = load_BulkTranscriptome(dataset, gene_num, qval).loc[ list(CCN.index), :]

# edf = load_BulkTranscriptome(dataset, gene_num=gene_num, qval=qval, NetworkSelection=False, standardize="StandardScaler")

# cellprop = pd.DataFrame( load_celltypeProportion(dataset, celltype=None, standardize="StandardScaler") ).loc[list(CCN.index), :]



## experiment
CCN_perf, CCN_result, CCN_weight = Leaveoneout_CV(CCN, list(clinical["response"]), feature_name="CCN", ML="LogisticRegression")
perf_list.append(CCN_perf)
result_list.append(CCN_result)

'''
bulk_perf, bulk_result, _ = Leaveoneout_CV(bulkseq_data, list(clinical["response"]), feature_name="bulkCommgenesNet", ML="LogisticRegression")
perf_list.append(bulk_perf)
result_list.append(bulk_result)

for name, control_geneset in zip(geneset_name, geneset):
    exist_geneset = [gene for gene in control_geneset["gene_id"].to_list() if gene in list(edf.columns)]
    control_edf = edf.loc[:, exist_geneset ]
    perf, result, _ = Leaveoneout_CV( control_edf, clinical["response"].to_list(), feature_name=name )

    perf_list.append( perf )
    result_list.append( result )

cellprop_perf, cellprop_result, _ = Leaveoneout_CV(cellprop, list(clinical["response"]), feature_name="cellprop", ML="LogisticRegression")
perf_list.append(cellprop_perf)
result_list.append(cellprop_result)
'''

## result 
perf = pd.concat(perf_list, axis=0)
print(perf)

result = pd.concat(result_list, axis=0)

perf.to_csv(result_dir + "/perf_" + dataset + ".txt", sep="\t", index_label="Feature")
result.to_csv(result_dir + "/result_" + dataset + ".txt", sep="\t", index_label="patient")
CCN_weight.to_csv(result_dir + "/CCNweight_" + dataset + ".txt", sep="\t", index_label="feature")



