import os
import shutil
import subprocess

from .JobDistributor import JobDistributor

class SGEJobDistributor(JobDistributor):
    '''Run jobs parallel on a sun grid engine.'''
    def __init__(self, data_set_name, script_name, script_arguments=[]):
        self.data_set_name = data_set_name
        self.script_name = script_name
        self.script_arguments = script_arguments
  
    def run(self, num_jobs, time='48:00:00', mem_free_GB=4, scratch_space_GB=1,
            architecture='linux-x64', hold_jid=None, keep_job_output_path=True):
        data_set_path = self.create_new_data_set(self.data_set_name)
        job_output_path = os.path.join(data_set_path, "job_outputs")
        
        if not keep_job_output_path and os.path.exists(job_output_path): # Clear the job output path
            shutil.rmtree(job_output_path)
 
        if not os.path.exists(job_output_path):
            os.mkdir(job_output_path)
 
        hj = ['-hold_jid', hold_jid] if hold_jid else []

        qsub_command = ['qsub',
                        '-cwd'] + hj \
                     + ['-N', self.script_name.split('/')[-1],
                        '-t', '1-{0}'.format(num_jobs),
                        '-l', 'h_rt={0}'.format(time),
                        '-l', 'mem_free={0}G'.format(mem_free_GB),
                        '-l', 'scratch={0}G'.format(scratch_space_GB),
                        '-l', 'arch=linux-x64',
                        '-o', job_output_path,
                        '-e', job_output_path,
                        './job_scripts/run_SGE_job.sh',
                        self.script_name,
                        data_set_path] \
                        + self.script_arguments \
                        + [num_jobs]
        
        subprocess.check_call(qsub_command)

    @staticmethod
    def get_sub_list(total_list, num_jobs, job_id):
        '''Get the sub-list responsible for a job to handle.'''
        return [total_list[i] for i in range(len(total_list)) if i % num_jobs == job_id]
