# NCBI Tool Kit

A tool kit for downloading and curating collections of genomes retrieved from the  National Center for Biotechology Information's public database, [GenBank](https://www.ncbi.nlm.nih.gov/).

   - Automatically synchronize your local collection with the [latest assembly versions](https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/#current).
   - Give FASTAs useful names based on information avaialable in the [assembly summary file](ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt) and the [taxonomy dump file](ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_readme.txt).

Requires [rsync](https://rsync.samba.org/).  Tested only with rsync version 3.1.2  protocol version 31.

## Installation

```
pip install ncbitk
ncbitk --help
```

Or simply clone this repository:

```
mkdir -p $HOME/projects/NCBITK && cd $HOME/projects/NCBITK
git clone https://github.com/andrewsanchez/NCBITK.git
pip install -r requirements.text
python run.py --help
```

## Usage

[![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg?style=flat-square)](http://makeapullrequest.com)
