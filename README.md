# PRS-Net
## About
This repository contains the code and resources of the following paper:

PRS-Net: Interpretable polygenic risk scores via geometric learning

## Overview of the framework
PRS-Net is an interpretable genomic deep learning-based approach designed to effectively model the nonlinearity of the biological system and deliver precise and interpretable polygenic risk scores.
<p align="center">
<img  src="framework.png"> 
</p>


## **Data preprocessing**
### **Setup environment**
Setup the required environment using `environment_data.yml` with Anaconda. While in the project directory run:

    conda env create -f environment_data.yml

Activate the environment

    conda activate PRS-Net_data
**Step 1: Arrange GWAS summary statistic data**

Ensure the columns in your GWAS summary statistic data are ordered as follows, separated by tabs (\t): CHR, BP, SNP, A1, A2, N, SE, P, OR, INFO, MAF, BETA, where A1 stands for the effect allele of the SNP and A2 stands for the non-effect allele of the SNP.

**Step 2: Run preprocessing script**

Execute the 1_data_preprocess.sh script with the following command:

    ./1_data_preprocess.sh your_base_data_path your_target_data_path your_phenotype_data_path your_ld_data_path your_output_path

* **your_base_data_path:** the path of the gwas summarty stastistics file derived from Step 1

* **your_target_data_path:** the path (without postfix) of the target data in PLINK format (.bim, .fam, .bed)

* **your_phenotype_data_path:** the path of the phenotype data separated by blanks (Sample_ID Sample_ID Phenotype_label)

* **your_ld_data_path:** the path (without postfix) of the data used for LD calculating in PLINK format

* **your_output_path:** the path of the output files

**Step 3: Run PRS-Net data generation script**

Execute the gwas_data_preprocess.sh script with the following command:

    ./2_data_generation.sh your_preprocessed_data_path your_phenotype_data_path your_output_path

* **your_preprocessed_data_path:** the path of the preprocessed files derived from Step 2

* **your_phenotype_data_path:** the path of the phenotype data separated by tabs (ID \t ID \t Label)

* **your_output_path:** the path of the output files




## **Train**
### **Setup environment**
Setup the required environment using `environment.yml` with Anaconda. While in the project directory run:

    conda env create -f environment.yml

Activate the environment

    conda activate PRS-Net


We upload an example dataset at [https://cloud.tsinghua.edu.cn/f/30a1c8df0c024fefa3c1/?dl=1](https://cloud.tsinghua.edu.cn/f/30a1c8df0c024fefa3c1/?dl=1).

Unzip the file and put it in example_dataset/ and excecute the train.py with the following command:

    python train.py --data_path ../example_dataset/ --dataset ad_eur


## License
PRS-Net is licensed under the Apache License, Version 2.0: http://www.apache.org/licenses/LICENSE-2.0.

