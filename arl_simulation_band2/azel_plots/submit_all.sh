#!/bin/bash
#!

cd case3
sbatch --dependency=singleton --job-name=IMAGING submit_p3.slurm
cd -

cd case5/100s
sbatch --dependency=singleton --job-name=IMAGING submit_p3.slurm
cd -

cd case5/30s
sbatch --dependency=singleton --job-name=IMAGING submit_p3.slurm
cd -

cd case5/10s
sbatch --dependency=singleton --job-name=IMAGING submit_p3.slurm
cd -

cd case6/m75
sbatch --dependency=singleton --job-name=IMAGING submit_p3.slurm
cd -

cd case6/0
sbatch --dependency=singleton --job-name=IMAGING submit_p3.slurm
cd -

cd case6/p15
sbatch --dependency=singleton --job-name=IMAGING submit_p3.slurm
cd -

cd case7
sbatch --dependency=singleton --job-name=IMAGING submit_p3.slurm
cd -

cd case10/delta-ra
sbatch --dependency=singleton --job-name=IMAGING submit_p3.slurm
cd -

cd case10/delta-dec
sbatch --dependency=singleton --job-name=surface submit_p3.slurm
cd -

cd case2
sbatch --dependency=singleton --job-name=IMAGING submit_p3.slurm
cd -


