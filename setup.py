#!/usr/bin/env python3
from setuptools import setup
import RMseq
import os


LONG_DESCRIPTION = 'Given a set of paired-end reads, a reference sequence ' +\
                   'and a reference protein (both in .fasta format), map ' +\
                   'the reads to the reference and determine the pool ' +\
                   'of mutated alleles in the short reads.'

if os.path.exists('README'):
    LONG_DESCRIPTION = open('README').read()

setup(
    name = 'rmseq',
    version = RMseq.__version__,
    description = RMseq.__description__,
    long_description=LONG_DESCRIPTION,
    classifiers = ['Development Status :: 3 - Alpha',
                   'License :: OSI Approved :: GNU Affero General ' +
                   'Public License v3 or later (AGPLv3+)',
                   'Programming Language :: Python :: 3.5',
                   'Topic :: Scientific/Engineering :: Bio-Informatics',
                   'Topic :: Scientific/Engineering :: Medical Science Apps.',
                   'Intended Audience :: Science/Research'],
    keywords = ['resistance',
                'mutation'],
    download_url = RMseq.__download_url__,
    author = RMseq.__author__,
    author_email = RMseq.__author_email__,
    license = RMseq.__license__,
    packages = ['RMseq'],
    scripts = ['RMseq/rmseq',
               'RMseq/amplicon-effect.py',
               'RMseq/RM-seq.pl'],
    include_package_data = True,
    install_requires = ['plumbum>=1.6.3',
                        'biopython>=1.69',
                        'numpy>=1.13.1']
    )
