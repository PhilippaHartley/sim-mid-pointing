#!/usr/bin/env bash

python ../pointing_simulation_distributed.py --context singlesource --rmax 1e5 --flux_limit 1.0 --ngroup 128 \
--static_pe 1.0 --dynamic_pe 0.0 --nworkers 8 --show True --seed 18051955