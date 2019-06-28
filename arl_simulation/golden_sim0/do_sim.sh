#!/usr/bin/env bash

python ../pointing_simulation_distributed.py --context s3sky --rmax 1e5 --flux_limit 0.3 --ngroup 8 \
--nworkers 8 --show True --seed 18051955  --pbtype MID_GAUSS --memory 16  --integration_time 10 \
--use_agg True --time_series wind --time_chunk 180 --reference_pointing True | tee pointing_simulation.log