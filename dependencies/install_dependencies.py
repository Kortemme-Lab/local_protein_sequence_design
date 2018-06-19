#!/usr/bin/env python2.7
'''
Install dependencies
'''

import os
import subprocess


def install_psipred():
    '''Install PSIPRED.'''
    subprocess.check_call(['wget', 'http://bioinfadmin.cs.ucl.ac.uk/downloads/psipred/psipred.4.01.tar.gz'])
    subprocess.check_call(['tar', '-xzf', 'psipred.4.01.tar.gz' ])
    
    print "PSIPRED is installed. Please set the runpsipred_single script properly."


if __name__ == '__main__':

    if not os.path.exists('dependencies'):
        os.mkdir('dependencies')

    os.chdir('dependencies')

    install_psipred()
