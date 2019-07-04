#!/bin/bash
#!
#! Dask job script for CSD3
#! Vlad Stolyarov
#!

#!#############################################################
#!#### Modify the options in this section as appropriate ######
#!#############################################################

#! sbatch directives begin here ###############################
#! Name of the job:
#SBATCH -J GLD0
#! Which project should be charged:
#SBATCH -A SKA-SDP-SL2-CPU
#! How many whole nodes should be allocated?
#SBATCH --nodes=4
#! How many (MPI) tasks will there be in total? (<= nodes*32)
#SBATCH --ntasks=32
#! Memory limit
#SBATCH --mem 128000
#! How much wallclock time will be required?
#SBATCH --time=6:00:00
#! What types of email messages do you wish to receive?
#SBATCH --mail-type=ALL
#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
##SBATCH --no-requeue
#! Rase the job priority for --time < 1h and --nodes <= 3
##SBATCH --qos=INTR
#! Do not change:
#SBATCH -p skylake
#! What types of email messages do you wish to receive?
#SBATCH --mail-type=FAIL
#! Ask for exclusive access
#SBATCH --exclusive
#! Same switch
#SBATCH --switches=1
#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
##SBATCH --no-requeue
#! Optionally modify the environment seen by the application
#! (note that SLURM reproduces the environment at submission irrespective of ~/.bashrc):
#module purge                               # Removes all modules still loaded

#! Set up python
. $HOME/rds/hpc-work/software/arlenv/bin/activate
#export PYTHONPATH=$PYTHONPATH:$ARL
#echo "PYTHONPATH is ${PYTHONPATH}"
echo -e "Running python: `which python`"
echo -e "Running dask-scheduler: `which dask-scheduler`"

module load slurm

JOBID=${SLURM_JOB_ID}
echo ${SLURM_JOB_NODELIST}

cd $SLURM_SUBMIT_DIR
echo -e "Changed directory to `pwd`.\n"

#! Create a hostfile:
scontrol show hostnames $SLURM_JOB_NODELIST > hostfile.$JOBID

hostIndex=0
for host in `cat hostfile.$JOBID`; do
    ibhost=${host}
    echo "Working on ${ibhost} ...."
    if [ "$hostIndex" = "0" ]; then
        scheduler=${ibhost}
        echo "run dask-scheduler"
        echo "ssh ${ibhost} dask-scheduler --host ${scheduler} --port=8786 &"
        ssh ${host} dask-scheduler --host ${scheduler} --port=8786 &
#	hostIndex="1"        
	sleep 5
    fi	
#    else
    echo "run dask-worker"
    ssh $host dask-worker \
	    --host ${ibhost} \
	    --nprocs 8 \
	    --nthreads 1  \
	    --memory-limit 0.75 \
	    --local-directory /tmp \
	    ${scheduler}:8786 &
	    sleep 1
    hostIndex="1"
#    fi
done
echo "Scheduler and workers now running"

#! We need to tell dask Client (inside python) where the scheduler is running
export ARL_DASK_SCHEDULER=${scheduler}:8786
echo "Scheduler is running at ${scheduler}"

#root=${ARL}/workflows/scripts/pipelines/performance/pipelines-timings
#root=/home/vs237/rds/hpc-work/software/algorithm-reference-library/workflows/notebooks
root=/home/vs237/software.ext/sim-mid-pointing/arl_simulation/simulation27
CMD="python ../pointing_simulation_distributed.py --context s3sky --rmax 1e5 --flux_limit 0.3 --ngroup 8 --nworkers 16 --show True --seed 18051955  --pbtype MID_GAUSS --memory 8  --integration_time 10 --use_agg True --time_series wind --time_chunk 1800 --reference_pointing True | tee pointing_simulation.log"

echo "About to execute $CMD"
eval $CMD

exit

