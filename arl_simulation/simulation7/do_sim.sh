#!/usr/bin/env bash

python ../pointing_simulation_distributed.py --context s3sky --rmax 1e5 --flux_limit 0.3 --ngroup 128 \
--static_pe 0.0 --dynamic_pe 1.0 --nworkers 8 --show False --seed 18051955  --pbtype MID_GAUSS --memory 32 \
--use_agg False