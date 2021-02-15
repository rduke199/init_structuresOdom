#!/bin/bash
#SBATCH -t 14-00:00:00                                #Time for the job to run
#SBATCH --job-name=wt_psi4           #Name of the job
#SBATCH -N 1                                            #Number of nodes required
#SBATCH -n 8                                            #Number of cores needed for the job
#SBATCH -p CAL48M192_L
#SBATCH --account=col_cmri235_uksr  #Name of account to run under

#SBATCH --mail-type ALL                         #Send email on start/end
#SBATCH --mail-user rdu230@uky.edu      #Where to send email
#SBATCH --error=SLURM_wt_psi4_%j.err           #Name of error file
#SBATCH --output=SLURM_wt_psi4_%j.out  #Name of output file
#SBATCH --array=1-20%20

echo "Job $SLURM_JOB_ID running on SLURM NODELIST: $SLURM_NODELIST "
echo "Job started at : $(date)"

module load ccs/conda/psi4-1.3.2+ecbda83
ulimit -s unlimited
ulimit -l unlimited
export OMP_NUM_THREADS=1

line=$(echo $SLURM_ARRAY_TASK_ID )
runfolder=$(sed -n "$line"p folders_to_run.txt)
echo $runfolder
cd $runfolder
time python ipfitting.py *xyz
