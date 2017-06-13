===============
 NCBI Tool Kit
===============

A tool kit for downloading and curating collections of genomes retrieved from the  National Center for Biotechology Information's public database,  `GenBank <https://www.ncbi.nlm.nih.gov/>`_.  NCBITK currenlty only supports downloading bacteria genomes.

* Automatically synchronize your local collection with the `latest assembly versions <https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/#current>`_.
* Give FASTAs useful names based on information avaialable in the `assembly summary file <ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt>`_ and the `taxonomy dump file <ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_readme.txt>`_.

Requires `rsync <https://rsync.samba.org/>`_.  Tested only with rsync version 3.1.2  protocol version 31.

==============
 Installation
==============

Using `pip <https://packaging.python.org/installing/>`_::

  pip install ncbitk

Or::

  git clone https://github.com/andrewsanchez/NCBITK.git
  python setup.py install

Regardless of which installation method you choose, I recommend using a `virtual environment <http://python-guide-pt-br.readthedocs.io/en/latest/dev/virtualenvs/>`_.

=============
 Basic Usage
=============

.. image:: https://github.com/andrewsanchez/NCBITK/raw/master/NCBITK/docs/NCBITK-Workflow.png

Download all GenBank bacteria::

  ncbitk [directory] --update

If you have already run NCBITK, the above will also update your local collection, i.e. remove old genomes no longer in the assembly summary and download the latest assembly versions.

Download only E. coli genomes::

  ncbitk [directory] --species Escherichia_coli --update

Note that in the above command, the list of strings given to the ``--species`` option must match exactly a species directory at ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/

Get the status of your collection::

  ncbitk [directory] --status

This will tell you how many genomes you have, what is missing from your collection, and how many deprecated genomes are present.


.. image:: https://img.shields.io/badge/PRs-welcome-brightgreen.svg?style=flat-square
