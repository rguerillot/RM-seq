from setuptools import setup
import RMseq
import os


LONG_DESCRIPTION = 'RM-seq is a bioinformatics tool for for ' +\
                   'assessing resistance mutations from PE short-reads.'


setup(
    name = 'rmseq',
    version = RMseq.__version__,
    description = RMseq.__description__,
    long_description=long_description,
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
    scripts = ['RMseq/rmseq'],
    include_package_data = True,
    install_requires = []
    )