#!/bin/bash

# run as sudo
# takes one argument for organism you want
# must match a directory name at ftp.ncbi.nlm.nih.gov::genomes/genbank/bacteria
# e.g. Escherichia_coli
# pass "\*" without quotes to sync with entire bacteria database

organism="$1"
ftplocation="bacteria/"$organism"/latest_assembly_versions/"
ftpmirror="./"$organism""
justfastas="$ftpmirror"_local

mkdir -p "$ftpmirror"
mkdir -p "$justfastas"

rsync -iPrLtm -f="- **unplaced_scaffolds**" -f="+ *.fna.gz" -f="+ */" -f="- *" ftp.ncbi.nlm.nih.gov::genomes/genbank/"$ftplocation" "$ftpmirror"

sudo find "$ftpmirror" -type f -exec cp -t "$justfastas" -- {} +

# sudo python /home/truthling/MGGen/NCBI_tools/rename.py

# --exclude-from=FILE     read exclude patterns from FILE
# --include-from=FILE     read include patterns from FILE
# -c, --checksum              skip based on checksum, not mod-time & size
# -i, --itemize-changes       output a change-summary for all updates
# -P, --partial, --progress   show progress and put partially downloaded files in a folder
# -r, --recursive             recurse into directories
# -L, --copy-links            transform symlink into referent file/dir
# -t, --times                 preserve modification times
# -m, --prune-empty-dirs      prune empty directory chains from file-list
# -n, --dry-run               perform a trial run with no changes made
