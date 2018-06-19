import os

class JobDistributor:
  '''Base class of job distributors.'''
  def __init__(self):
    pass

  def create_new_data_set(self, name, data_dir='data'):
    # Create the directory to save all data
    
    if not os.path.exists(data_dir): os.makedirs(data_dir)
  
    # Create a new data set directory with the name
  
    data_set_path = os.path.join(data_dir, name)
    if not os.path.exists(data_set_path) : os.makedirs( data_set_path )
  
    return data_set_path
