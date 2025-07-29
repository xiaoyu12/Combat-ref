#!/bin/bash
Rscript sim_DEpipe.R 1 1 12 &
Rscript sim_DEpipe.R 1 2 12 &
Rscript sim_DEpipe.R 1 3 12 &
Rscript sim_DEpipe.R 1 4 12 &
Rscript sim_DEpipe.R 1.5 1 12 &
Rscript sim_DEpipe.R 1.5 2 12 &

Rscript sim_DEpipe.R 1.5 3 12 &
Rscript sim_DEpipe.R 1.5 4 12 &
Rscript sim_DEpipe.R 2 1 12 &
Rscript sim_DEpipe.R 2 2 12 &
Rscript sim_DEpipe.R 2 3 12 &
Rscript sim_DEpipe.R 2 4 12 &
Rscript sim_DEpipe.R 2.4 1 12 &
Rscript sim_DEpipe.R 2.4 2 12 &
Rscript sim_DEpipe.R 2.4 3 12 &
Rscript sim_DEpipe.R 2.4 4 12 &


wait
echo "All simulations complete."
