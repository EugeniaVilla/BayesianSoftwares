#PBS -l nodes=1:ppn=6,walltime=24:00:00 -q gigat

# Set the job name
#PBS -N synthetic_experiments_gutierrez

# Set the output file and merge it to the sterr
#PBS -o out-hostname-XyZ-N1x1-qsub.txt
#PBS -j oe
#PBS -e out-hostname-XyZ-N1x1.txt

# Send emails
#PBS -M mario.beraha@polimi.it
#PBS -m abe


#------------------SETTING THE ENVIRONMENT----------------------------#
# Start the job in the current directory (PBS starts in the home folder)
cd ${PBS_O_WORKDIR}

# export MY_ENVIRONMENTAL_VARIABLE=value

export mkPrefix=/u/sw
source $mkPrefix/etc/profile
module load gcc-glibc/6
module load R


#-------------------RUN THE EXECUTABLE---------------------------------#

Rscript four_groups.R &> four_groups_log.txt