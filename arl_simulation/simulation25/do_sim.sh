#!/usr/bin/env bash
cp ../pointing_simulation.py .
python pointing_simulation.py --context singlesource --rmax 1e5 --flux_limit 1.0 --ngroup 128 \
--static_pe 0.0 --dynamic_pe 1.0 --nworkers 8 --show True --seed 18051955 --scale 1.0 1.0 --use_agg True \
--use_radec False --integration_time 10 --time_range -6 6 --pbtype MID_GAUSS \
--memory 32 --time_series wind
