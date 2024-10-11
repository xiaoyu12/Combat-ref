# Combat-ref
We introduce ComBat-ref, a refined method of batch effect correction that enhances the statistical power and reliability of differential expression analysis in RNA-seq data. Building on the foundations of ComBat-seq, ComBat-ref employs a negative binomial model to adjust count data but innovates by using a pooled dispersion parameter for entire batches and preserving count data for the reference batch. Our method demonstrated superior performance in both simulated environments and real datasets. By effectively mitigating batch effects while maintaining high detection power, ComBat-ref proves to be a robust tool for enhancing the accuracy and interpretability of RNA-seq data analyses.

## Run the tests
Run the simulation test with the following command, where batch_FC = 2, disp_FC = 4, and a total of 12 samples (3 for each batch and biological condition). You can change these parameters to compare batch effect correction results by Combat-ref and other methods. 
```
Rscript simulations/sim_DEpipe.R 2 4 12
```

To run the test for GFRN data, use the following script.
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
