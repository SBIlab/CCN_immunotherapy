import pandas as pd
import numpy as np
import os
import sklearn.metrics as metrics
from sklearn.preprocessing import StandardScaler
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests


def nonres_res_partition(dataset):
    clinical_path = os.getcwd() + "/../../data/CCN_construction/1_CIBERSORTx_input/" + dataset + "/"

    if dataset == "Liu":
       clinical = pd.read_csv(filepath_or_buffer=clinical_path + "clinical.txt", sep="\t", index_col="ID")
       clinical = clinical.loc[ clinical["response"] != "MR", :]
       clinical["respond"] = [1 if response == "PR" or response == "CR" else 0 for response in clinical["response"]]
       clinical = clinical.loc[:, ["respond"]].rename(columns={"respond": "response"})

    elif dataset == "Gide":
        clinical = pd.read_csv(filepath_or_buffer=clinical_path + "clinical.txt", sep="\t", index_col="ID")
        clinical = clinical.loc[ : , ["flag"]].rename(columns={"flag": "response"})
    
    elif dataset == "Mariathasan":
        clinical = pd.read_csv(filepath_or_buffer=clinical_path + "pData.txt", index_col="sample_id")
        clinical = clinical.loc[ clinical["Best Confirmed Overall Response"] != "NE", ["Best Confirmed Overall Response"]]
        clinical["response"] = [1 if response == "PR" or response == "CR" else 0 for response in clinical["Best Confirmed Overall Response"]]
        clinical = clinical.loc[:, ["response"]]
    
    elif dataset == "Kim":
        clinical = pd.read_csv(filepath_or_buffer=clinical_path + "clinical.txt", sep="\t", index_col=0)
        clinical = clinical.loc[:, ["response"]]
    
    elif dataset == "Hugo":
        clinical = pd.read_csv(filepath_or_buffer=clinical_path + "clinical.txt", sep="\t", index_col="geo_accession")
        clinical = clinical.loc[clinical["biopsy_time"] == "pre-treatment",:]
        clinical["response"] = [0 if response == "Progressive Disease" else 1 for response in clinical["respond"]]
        clinical = clinical.loc[ :, ["response"] ]

    elif dataset == "VanAllen":
        clinical = pd.read_csv(filepath_or_buffer=clinical_path + "clinical.txt", sep="\t", index_col=0)
        clinical = clinical.loc[:, ["response"]]

    elif dataset == "Cho":
        clinical = pd.read_csv(filepath_or_buffer=clinical_path + "clinical.txt", sep="\t", index_col="geo_accession")
        clinical = clinical.loc[:,["Best_response"]].rename(columns={"Best_response": "respond"})
        clinical["response"] = [ 0 if response == "PD" or response == "SD" else 1 for response in clinical["respond"] ]
        clinical = clinical.loc[:, ["response"]]
    
    elif dataset == "Jung":
        clinical = pd.read_csv(filepath_or_buffer=clinical_path + "clinical.txt", sep="\t", index_col="geo_accession")
        clinical = clinical.loc[:, ["response"]]    
    
    elif dataset == "PratMelanoma":
        clinical = pd.read_csv(filepath_or_buffer=clinical_path + "clinical.txt", sep="\t", index_col="samples")
        clinical = clinical.loc[:, ["best.response"]]
        clinical["response"] = [1 if response == "PR" or response == "CR" else 0 for response in clinical["best.response"]]
        clinical = clinical.loc[:, ["response"]]


    return clinical



def NetworkPropagated_enriched_CommPathway(dataset, gene_number, cutoff):
    if dataset == "Liu": 
        drug_target = "PD1"
    elif dataset == "Gide": 
        drug_target = "PD1_CTLA4"
    elif dataset == "Mariathasan": 
        drug_target = "PD-L1"
    elif dataset == "Kim":
        drug_target = "PD1"
    elif dataset == "Hugo":
        drug_target = "PD1"
    elif dataset == "Cho":
        drug_target = "PD1"
    elif dataset == "Jung":
        drug_target = "PD1_PD-L1"
    elif dataset == "PratMelanoma": 
        drug_target = "PD1"
    elif dataset == "VanAllen": 
        drug_target = "CTLA4"

   
    dataset_edf = pd.read_csv(os.getcwd() + "/../../data/ML/BulkTranscriptome/" + dataset + "/nonresponder_BulkTranscriptome.txt",
                             sep="\t")
    dataset_genelist = list(dataset_edf.columns)
    

    propa_score = pd.read_csv(os.getcwd() + "/../../data/ML/NetworkPropagatedScores/" + drug_target + ".txt", sep="\t")
    propa_score = propa_score.sort_values(by=["propagate_score"], ascending=False)
    propa_score = propa_score.loc[ propa_score["gene_id"].isin(dataset_genelist), :]
    high_propascored_gene = propa_score["gene_id"].to_list()[0:gene_number]

    CommPathway_genelist = pd.read_csv(os.getcwd() + "/../../data/ML/Communication_genelist/pathway/CommunicationPathway_genelist.txt", sep="\t")
    Commpathway = {}
    for key, value in zip(CommPathway_genelist["pathway"].values, CommPathway_genelist["gene_id"].values): 
        Commpathway[key] = value.split(",")
    
    geom_pvalue = []
    for key in list(Commpathway.keys()):
        M = len(propa_score["gene_id"])
        n = len(high_propascored_gene)
        N = len(Commpathway[key])
        k = len(set(high_propascored_gene) & set(Commpathway[key]))
        geom_pvalue.append( hypergeom.sf(k-1, M, n, N) )
    _, adj_pvalue, _, _ = list( multipletests(geom_pvalue, method="fdr_bh") )
    Commpathway_df = pd.DataFrame([list(Commpathway.keys()), geom_pvalue, adj_pvalue], index=["pathway", "pvlaue", "adj_pvalue"]).T
    enriched_Commpathway_df = Commpathway_df.loc[Commpathway_df["adj_pvalue"] < cutoff, :]

    return enriched_Commpathway_df.sort_values(by=["adj_pvalue"])



def load_celltypeProportion(dataset, celltype=None, standardize="StandardScaler"):
    path = os.getcwd() + "/../../data/CCN_construction/2_CIBERSORTx_output/" + dataset + "/CIBERSORTxGEP_NA_Fractions-Adjusted.txt"
    cellprop = pd.read_csv(path, index_col=0, sep="\t").drop(columns=["P-value", "Correlation", "RMSE"])

    if standardize == "StandardScaler":
        cellprop = pd.DataFrame( StandardScaler().fit_transform(cellprop), columns=list(cellprop.columns), index=list(cellprop.index) )
        
    return cellprop



def load_CCN(dataset, gene_num, p_cut_off, standardize="StandardScaler", NetworkSelection=True):
    path = os.getcwd() + "/../../data/CCN_construction/5_CCN/" + dataset + "/"
    
    
    CCN_list = [pd.read_csv(path + "nonresponder_CCN.txt", sep="\t", index_col=0), 
                pd.read_csv(path + "responder_CCN.txt", sep="\t", index_col=0) ]
    CCN = pd.concat(CCN_list, axis=0)
    

    if standardize == "StandardScaler":
        CCN = pd.DataFrame( StandardScaler().fit_transform(CCN), columns=list(CCN.columns), index=list(CCN.index) )
    

    if NetworkSelection:
        enriched_Commpathway_df = NetworkPropagated_enriched_CommPathway(dataset, gene_number=gene_num, cutoff=p_cut_off)
        column_list = [ col for col in CCN.columns if col.split(".")[0] in list(enriched_Commpathway_df["pathway"]) ]
        CCN = CCN.loc[:, column_list]
    
    print(dataset + " CCN loaded: ")

    return CCN





