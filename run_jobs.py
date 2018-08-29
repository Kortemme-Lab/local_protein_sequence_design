#!/usr/bin/env python3

'''
Run jobs for binding site sequence design. Before launching this 
wrapper script. you should write your costumized script the specify what
to do.
Then run your scripts with this wrapper.

Usage:
    ./run_jobs.py <name> <job-script> [options]

Arguments:
    <name>
        Name of the job.

    <job-script>
        Job script for the analysis.

Options:
    --job-distributor=<JD>, -d=<JD>  [default: sequential]
        Job distributor that runs the jobs.

    --hold-jid=<JID>, -h=<JID>  
        For SGEJobDistributor, hold the job until some jobs are finished.

    --keep-job-output-path, -k
        For SGEJobDistributor, do not clear the job output path.

    --job-script-arguments=<JA>, -a=<JA>  [default: ]
        Arguments passed to the job script. The job script will run
        as:
            job-script data-set-path job-script-arguments

    --num-jobs=<NJ>, -n=<NJ>
        Number of jobs for parallel run.
'''

import docopt

import local_protein_sequence_design as LPSD


if __name__ == '__main__':

    arguments = docopt.docopt(__doc__)
    
    # Convert job-script arguments into a list
    
    job_script_arguments = arguments['--job-script-arguments'].split()
    
    # Initialize the job distributor and run the job
    
    job_distributor = None
    
    if arguments['--job-distributor'] == 'sequential':
        job_distributor = LPSD.job_distributors.SequentialJobDistributor(arguments['<name>'],
                          arguments['<job-script>'], job_script_arguments)
        
        job_distributor.run()
    
    elif arguments['--job-distributor'] == 'SGE':
        job_distributor = LPSD.job_distributors.SGEJobDistributor(arguments['<name>'],
                          arguments['<job-script>'], job_script_arguments)
    
        num_jobs = arguments['--num-jobs'] if arguments['--num-jobs'] else 1
    
        job_distributor.run(num_jobs, hold_jid=arguments['--hold-jid'], keep_job_output_path=arguments['--keep-job-output-path'])
    
    else:
        raise IOError('Unknown job distributor: {0}'.format(arguments['--job-distributor']))
