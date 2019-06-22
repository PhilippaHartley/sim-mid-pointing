#!/usr/bin/env bash

python ../pointing_simulation_distributed.py --context singlesource --rmax 1e5 --flux_limit 1.0 --ngroup 8 \
--static_pe 0.0 --dynamic_pe 1.0 --nworkers 16 --show True --seed 18051955 --time_chunk 1800 --use_agg True \
| tee pointing_simulation.log