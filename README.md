# CCN_immunotherapy
## Description
This is source code for constructing cell-cell communication networks using patient's bulk tumor transcriptome and prediction of patient response to immune checkpoint inhibitors (ICIs) using the patient's cell-cell communication network (CCN)

## Packages required
1. construction of CCN
- R (v 4.0.2)
- preprocessCore (v 1.52.1)
- CellChat (v 1.1.0)
- dplyr (v 1.0.10)
- Matrix (v 1.5.1)
- reticulate (v. 1.26)

2. response prediction using CCN
- python (v 3.6.13)
- pandas (v 1.1.5)
- numpy (v 1.19.2)
- sklearn (v. 0.24.2)
- scipy (v 1.5.2)

## Installation

    git clone https://github.com/SBIlab/CCN_immunotherapy


## Usage with examples
### 1. Construction of cell-cell communication networks
- **pre-requirement: cell type-specific gene expression profiles should be prepared before running and store under ./data/CCN_construction/2_CIBERSORTx_output/"dataset"/
- input: cell type-specific gene expression profiles
- To construct CCN of each patient, run R codes deposited in ./code/CCN_construction (setting working directory and dataset variable at the beginning of each code)
- CCN will be stored under ./data/CCN_construction/5_CCN/"dataset"/
- running examples: 
  
      cd code/CCN_construction
      
      # preprocessing cell type-specific gene expression profiles
      Rscript 2to3_Processing_CIBERSORTxOutput_todo_CellChat.R    # results will be stored in 3_processed_celltype_gene_expression

      # calculate communication strength (i.e. communication probabilities)
      Rscript 3to4_Calculate_commprob_byCellChat.R        # results will be stored in 4_Commprob_calculated_CellChatObj

      # generate patient-specific cell-cell communication networks
      Rscript 4to5_CCN_extract_from_CellChatObj.R         # results will be stored in 5_CCN

### 2. Response prediction using CCN
- input: communication strength between cells of each patient (i.e. CCN)
- To make leave-one-out cross-validation using CCN, run ./code/ML/main.py under the ./code/ML/ directory
- Results will be stored at ./result/ML/"dataset"/
- running examples

      python main.py

## Example dataset
Jung et al., 2019, Nat. Commun. 10, 4278.
