#!/usr/bin/env bash

python ../pointing_simulation_distributed.py --context s3sky --rmax 1e5 --flux_limit 3 --ngroup 8 \
--nworkers 16 --show True --seed 18051955  --pbtype MID_GAUSS --memory 8  --integration_time 100 \
--use_agg True --time_series wind --time_chunk 1800 --reference_pointing True | tee pointing_simulation.log