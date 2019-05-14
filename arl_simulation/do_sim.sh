#!/usr/bin/env bash
python pointing_simulation.py --context singlesource --rmax 1e5 --flux_limit 1.0 --ngroup 128 --static_pe 0.0 --dynamic_pe 1.0 --nworkers 16
python pointing_simulation.py --context singlesource --rmax 1e5 --flux_limit 1.0 --ngroup 128 --static_pe 1.0 --dynamic_pe 0.0 --nworkers 16
python pointing_simulation.py --context s3sky --rmax 1e5 --flux_limit 0.1 --ngroup 128 --static_pe 0.0 --dynamic_pe 1.0 --nworkers 16
python pointing_simulation.py --context s3sky --rmax 1e5 --flux_limit 0.1 --ngroup 128 --static_pe 1.0 --dynamic_pe 0.0 --nworkers 16
