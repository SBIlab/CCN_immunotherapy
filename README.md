# CCN_immunotherapy
## Description
This is source code for constructing cell-cell communication networks using patient's bulk tumor transcriptome and prediction of patient response to immune checkpoint inhibitors (ICIs) using the patient's cell-cell communication network (CCN)

## Requirements
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

## Usage
1. Construction of cell-cell communication networks
- input: cell type-specific gene expression profiles (store under ./data/CCN_construction/2_CIBERSORTx_output/"dataset"/)
- To construct CCN of each patient, run R codes deposited in ./code/CCN_construction (setting working directory at the beginning of each code)
- CCN will be stored under ./data/CCN_construction/5_CCN/"dataset"/

2. Response prediction using CCN
- input: communication strength between cells of each patient (i.e. CCN)
- To make leave-one-out cross-validation using CCN, run ./code/ML/main.py under the ./code/ML/ directory
- Results will be stored at ./result/ML/"dataset"/
