
generate_submission() {
    local outputfile="$1/submission.sh"
    local solver_path="$2"
    local input_file="$3"
    local job_name="$4"
    local ntasks="$5"
    local walltime=""
    local email=""

    if [ -z "$6" ]; then
        walltime="24:00:00"
    else
        walltime="$6"
    fi

    if [ -z "$7" ]; then
        echo "Error: Email address is required for job notifications."
        exit 1
    else
        email="$7"
    fi

    cat <<EOF > $outputfile
#!/bin/bash
# Job Name and Files (also --job-name)
#SBATCH -J $job_name

#Output and error (also --output, --error)
#SBATCH -o ./%x_%j.out
#SBATCH -e ./%x_%j.err

#Notification and type
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=$email

# Wall clock limit:
#SBATCH --time=$walltime
#SBATCH --no-requeue

#SBATCH --exclusive

# set number of tasks - should be an integer multiple of physical cores
# per node

#SBATCH --ntasks=$ntasks

echo "Solver: $solver_path"
echo "Input File: $input_file"

source $HOME_KADATH/.hpc_ci_runs/env_template.sh

mpirun ${solver_path}/bin/Release/solve $input_file .
EOF

}