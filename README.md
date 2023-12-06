# CCN_immunotherapy
## Description
This is source code for constructing cell-cell communication networks using patient's bulk tumor transcriptome and prediction of patient response to immune checkpoint inhibitors (ICIs) using the patient's cell-cell communication network (CCN)

## Packages required
We recommend users to utilize conda environment to easily manage some issues caused by software version. (refer to: https://www.anaconda.com/download)
R (v 4.0.2) and python (v 3.6.13) were used to develop for constructing cell-cell communication networks (CCN) and prediction of patient response to ICIs using CCN.
The working environment can be set up as follows (if conda environments was utilized):

      (base) conda create -n CCN python=3.6.13
      (base) conda activate CCN
      (CCN) conda install -c conda-forge r-base=4.0.2

For overall processes of constructing cell-cell communication networks (CCN) and response prediction using CCN, both R and python packages, described below, are required.
Before running codes, below packages should be installed.

### R
- R (v 4.0.2)
- dplyr (v 1.0.10)
- preprocessCore (v 1.52.1)
- Matrix (v 1.5.1)
- future (v 1.33.0) (or CellChat v 1.1.0)
- pbapply (v 1.7.2) (or CellChat v 1.1.0)
- reticulate (v. 1.26)


For installing R packages in R console,

      install.packages("dplyr", repos = "https://cloud.r-project.org")
      install.packages("BiocManager", repos = "https://cloud.r-project.org")    ## to install bioconductor package "preprocessCore"
      BiocManager::install("preprocessCore") 
      install.packages("Matrix", repos = "https://cloud.r-project.org")
      install.packages("future", repos = "https://cloud.r-project.org")      ## if you installed CellChat v 1.1.0, no need to install this package
      install.packages("pbapply", repos = "https://cloud.r-project.org")     ## if you installed CellChat v 1.1.0, no need to install this package
      install.packages("reticulate", repos = "https://cloud.r-project.org") 

if you encounter below error during running "Rscript 2to3_Processing_CIBERSORTxOutput_todo_CellChat.R"
    
    Error in normalize.quantiles(as.matrix(exp_df), copy = FALSE) :  
    ERROR; return code from pthread_create() is 22 
    Calls: %>% -> as.data.frame -> normalize.quantiles 
    Execution halted 
Then, install package "preprocessCore" again, as shown below

    BiocManager::install("preprocessCore", configure.args = c(preprocessCore = "--disable-threading"), force= TRUE, update=TRUE, type = "source")

### python
- python (v 3.6.13)
- pandas (v 1.1.5)
- numpy (v 1.19.2)
- sklearn (v. 0.24.2)
- scipy (v 1.5.2)
- statsmodels (v. 0.12.2)


For installing python modules using Anaconda

      conda install pandas=1.1.5     ## numpy is also installed during this step
      conda install scikit-learn=0.24.2    ## scipy is also installed during this step
      conda install statsmodels=0.12.2
Or using pip

    pip install pandas==1.1.5        ## numpy is also installed during this step
    pip install scikit-learn==0.24.2    ## scipy is also installed during this step
    pip install statsmodels==0.12.2 
      

## Installation

    git clone https://github.com/SBIlab/CCN_immunotherapy


## Usage with examples
### 1. Construction of cell-cell communication networks
- **pre-requirement: cell type-specific gene expression profiles should be prepared before running these codes and store under ./data/CCN_construction/2_CIBERSORTx_output/"dataset"/
- input: cell type-specific gene expression profiles
- To construct CCN of each patient, run R codes deposited in ./code/CCN_construction (setting dataset variable at the beginning of each code)
- CCN will be stored under ./data/CCN_construction/5_CCN/"dataset"/
- running example:
  
      cd CCN_immunotherapy/code/CCN_construction
      
      # preprocessing cell type-specific gene expression profiles
      Rscript 2to3_Processing_CIBERSORTxOutput_todo_CellChat.R    # results will be stored in 3_processed_celltype_gene_expression

      # calculate communication strength (i.e. communication probabilities)
      Rscript 3to4_Calculate_commprob_byCellChat.R        # results will be stored in 4_Commprob_calculated_CellChatObj

      # generate patient-specific cell-cell communication networks
      Rscript 4to5_CCN_extract_from_CellChatObj.R         # results will be stored in 5_CCN
- output directory and files:

    - data/CCN_construction/3_processed_celltype_gene_expression/"dataset"/ -> pre-processed cell type-specific gene expression profiles
    - data/CCN_construction/4_Commprob_calculated_CellChatObj/"dataset"/ -> CellChat object of each sample containing communication probabilities. Saved as "sampleName"_CellChatobj.RData
    - data/CCN_construction/5_CCN/"dataset"/ -> cell-cell communication networks of each sample, divided into responders (responder_CCN.txt) and non-responders (nonresponder_CCN.txt). 


### 2. Response prediction using CCN
- input: communication strength between cells of each patient (i.e. CCN)
- To make leave-one-out cross-validation using CCN, run ./code/ML/main.py under the ./code/ML/ directory
- Results will be stored at ./result/ML/"dataset"/
- running example

      cd CCN_immunotherapy/code/ML
      python main.py
- output directory and files:

    - result/ML/"dataset"/perf_"dataset".txt -> predictive performances of the CCN-based ML model
    - result/ML/"dataset"/result_"dataset".txt -> predictive probabilities of each sample by CCN-based ML model
    - result/ML/"dataset"/CCNweight_"dataset".txt -> feature weights

## Example dataset
Jung et al., 2019, Nat. Commun. 10, 4278.
