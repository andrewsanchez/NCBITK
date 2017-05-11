from setuptools import setup, find_packages
from os import path

setup(
    name = 'NCBITK',
    packages = ['NCBITK'], # this must be the same as the name above
    version = '1.0a1',
    license = 'MIT',

    here = path.abspath(path.dirname(__file__))
    with open(path.join(here, 'README.org'), encoding='utf-8') as f:
        long_description = f.read()

    description = "A tool kit for NCBI's GenBank",
    long_description=long_description,
    author = 'Andrew Sanchez',
    author_email = 'inbox.asanchez@gmail.com',
    url = 'https://github.com/andrewsanchez/NCBITK',
    keywords = 'NCBI bioinformatics',
    classifiers =[
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3.6',
        'Operating System :: POSIX :: Linux',
        'License :: OSI Approved :: MIT License',
        'Intended Audience :: Science/Research'
        'Environment :: Console',
        'Development Status :: 4 - Beta'
    ],
)
