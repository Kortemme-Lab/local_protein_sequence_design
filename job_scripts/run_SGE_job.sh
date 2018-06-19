#!/bin/bash
#
#$ -S /bin/bash
#$ -cwd

# This Bash script sets up the virtual environment and 
# runs python job scripts. Usage:
# ./run_job.sh job_script_command


# Set up the virtual environment
source venv/bin/activate

# Run the job_script_command
$@ ${SGE_TASK_ID}
