#-----------------SETTING THE REQUEST FOR THE HARDWARE ALLOCATION----#

#PBS -l nodes=1:ppn=1,walltime=24:00:00 -q gigat

# Set the job name
#PBS -N glaucoma_dp

# Set the output file and merge it to the sterr
#PBS -o out-hostname-XyZ-N1x1-qsub.txt
#PBS -j oe
#PBS -e out-hostname-XyZ-N1x1.txt


#------------------SETTING THE ENVIRONMENT----------------------------#
# Start the job in the current directory (PBS starts in the home folder)
cd /u/archive/dott/beraha/cmdstan-2.18.1

# export MY_ENVIRONMENTAL_VARIABLE=value

export mkPrefix=/u/sw
source $mkPrefix/etc/profile
module load gcc-glibc/7

module load eigen

# --------------- COMPILE THE MODEL -----------------------------#
make /u/beraha/time_dep_bnp/stan/glaucoma_dp &> /u/beraha/time_dep_bnp/stan/glaucoma_dp_compile_log.txt

#-------------------RUN THE EXECUTABLE---------------------------------#


/u/beraha/time_dep_bnp/stan/glaucoma_dp  \
  sample num_samples=2000 num_warmup=2000 thin=1 \
  data file=/u/archive/dott/beraha/data/time_dependent_bnp/data_imputed_regression_all_scaled_cmdstan.R \
  output file=/u/archive/dott/beraha/data/time_dependent_bnp/chain_dp_stan.csv  &> /u/beraha/time_dep_bnp/stan/glaucoma_dp_run_log.txt
