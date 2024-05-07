# Combat-ref

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
