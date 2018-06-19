import os
import json

binding_site_sequence_design_home = '.'

def load_site_settings(setting_file=None):
    '''Load the site settings into a dictionary'''
    if setting_file is None:
        setting_file = os.path.join(binding_site_sequence_design_home, 'site_settings/site_settings.json')
    
    with open(setting_file, 'r') as f:
        return json.load(f)
