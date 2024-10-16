# Highly Effective Batch Effect Correction Method for RNA-seq Count Data
We introduce ComBat-ref, a new method of batch effect correction that enhances the statistical power and reliability of differential expression analysis in RNA-seq data. Building on the foundations of ComBat-seq, ComBat-ref employs a negative binomial model to adjust count data but innovates by using a pooled dispersion parameter for entire batches and preserving count data for the reference batch. Our method demonstrated superior performance in both simulated environments and real datasets. By effectively mitigating batch effects while maintaining high detection power, ComBat-ref proves to be a robust tool for enhancing the accuracy and interpretability of RNA-seq data analyses.

![tpr_fpr_fc2_disp3](https://github.com/user-attachments/assets/7c0c2ecf-f195-4150-a7be-797bc5c8df7b)

![pca_variation_aggr](https://github.com/user-attachments/assets/96b2a30e-7ece-4186-bffc-005d9f7e700b)
## Run the tests
Run the simulation test with the following command, where batch_FC = 2, disp_FC = 4, and a total of 12 samples (3 replicates for each batch and biological condition). You can change these parameters to compare batch effect correction results by Combat-ref and other methods. In this test, the batch fold change (batch_FC) represents the mean batch effects, i.e., the ratio of the mean of batch 2 over that of batch 1; the dispersion batch effect (disp_FC) is the ratio of the dispersion of batch 2 over that of batch 1. The larger the batch_FC and disp_FC values, the more difficult it is to correct the batch effects. In the simulations, we assume the fold changes of biological signals to be 2.4. Assuming two biological conditions and two batches, 12 samples mean 3 replicates for each batch and biological condition.
```
Rscript simulations/sim_DEpipe.R 2 4 12
```
Visualize the simulation results using the following command:
```
Rscript simulations/plot_sim.R 2 4
```
It creates two plots comparing the true positive rates (TPR) and false positive rates (FPR) of different batch correction methods at p-value = 0.05 and FDR = 0.1.

To run the test for GFRN data, you can just use the following script.
```
Rscript real_data_application/gfrn_application.R
```
and
```
Rscript real_data_application/gfrn_DE.R
```

To run the test for NASA GenLab data, use the following script.
```
Rscript nasa_data/nasa_application.R
```
## Batch Correct Count Data
Run the "batch_correct.R" to correct the batch effects in RNA-seq count data. For example:
```
Rscript batch_correct.R gfrn.mat gfrn_samples.csv
```
The script takes two command-line arguments: the first is the count matrix and the second the sample description file. Please see the example "gfrn.mat" and "gfrn_samples.csv" for the file format. The corrected data is written to the "gfrn_corrected.mat" file.
## Cite
Highly Effective Batch Effect Correction Method for RNA-seq Count Data

Xiaoyu Zhang

bioRxiv 2024.05.02.592266; doi: https://doi.org/10.1101/2024.05.02.592266
