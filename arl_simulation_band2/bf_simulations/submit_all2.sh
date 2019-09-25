#!/bin/bash
#!

cd case5/10s
sbatch --dependency=singleton --job-name=TYPE1 submit_p3.slurm
cd -

cd case6/m75
sbatch --dependency=singleton --job-name=TYPE1 submit_p3.slurm
cd -

cd case7
sbatch --dependency=singleton --job-name=TYPE1 submit_p3.slurm
cd -

cd case10/delta-ra
sbatch --dependency=singleton --job-name=TYPE1 submit_p3.slurm
cd -

cd case10/delta-dec
sbatch --dependency=singleton --job-name=TYPE1 submit_p3.slurm
cd -

cd case2
sbatch --dependency=singleton --job-name=TYPE1 submit_p3.slurm
cd -


