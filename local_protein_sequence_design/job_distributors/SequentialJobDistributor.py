import subprocess

from .JobDistributor import JobDistributor


class SequentialJobDistributor(JobDistributor):
    '''Run jobs sequentially.'''
    def __init__(self, data_set_name, script_name, script_arguments=[]):
        self.data_set_name = data_set_name
        self.script_name = script_name
        self.script_arguments = script_arguments
  
    def run(self):
        data_set_path = self.create_new_data_set(self.data_set_name)
        subprocess.check_call([self.script_name, data_set_path] + self.script_arguments) 
