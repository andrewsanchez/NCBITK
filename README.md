# NCBI Tool Kit

A tool kit for downloading and curating collections of genomes retrieved from the  National Center for Biotechology Information's public database, [GenBank](https://www.ncbi.nlm.nih.gov/).

   - Automatically synchronize your local collection with the [latest assembly versions](https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/#current).
   - Give FASTAs useful names based on information avaialable in the [assembly summary file](ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt) and the [taxonomy dump file](ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_readme.txt).

Requires [rsync](https://rsync.samba.org/).  Tested only with rsync version 3.1.2  protocol version 31.

## Installation

Using [pip](https://packaging.python.org/installing/):

```
pip install ncbitk
```

Or:

```
git clone https://github.com/andrewsanchez/NCBITK.git
python setup.py install
```

Regardless of which installation method you choose, I recommend using a [virtual environment](http://python-guide-pt-br.readthedocs.io/en/latest/dev/virtualenvs/).

## Basic Usage

![NCBITK Workflow](/images/NCBITK-Workflow.png)

### Download all GenBank bacteria

```
ncbitk local-directory --update
```

If you have already run NCBITK, the above will also update your local collection, i.e. remove old genomes no longer in the assembly summary and download the latest assembly versions.

### Get the status of your collection

```
ncbitk local-directory --status
```

This will tell you how many genomes you have, what is missing from your collection, and how many deprecated genomes are present.


[![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg?style=flat-square)](http://makeapullrequest.com)
