#!/usr/bin/env python3

from distutils.core import setup

setup(
    name='local_protein_sequence_design',
    version='0.0.0',
    author='Xingjie Pan',
    author_email='xingjiepan@gmail.com',
    url='https://github.com/xingjiepan/local_protein_sequence_design',
    packages=[
        'local_protein_sequence_design',
    ],
    install_requires=[
        'numpy',
        'matplotlib',
        'docopt',
    ],
    extras_require = {
        'weblogo':  ['weblogo'],
        'pandas':  ['pandas'],
    },
    entry_points={
        'console_scripts': [
        ],
    },
    description='Design a local part of a protein',
    long_description=open('README.md').read(),
    classifiers=[
        'Programming Language :: Python :: 3',
        'Intended Audience :: Science/Research',
    ],
)
