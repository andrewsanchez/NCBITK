#!/bin/bash

# run as sudo

organism="Escherichia_coli" # set this to just "*" for all bacteria
location="bacteria/"$organism"/latest_assembly_versions/"
directory="/home/truthling/MGGen/"$organism"" # must be relative; will be created if doesn't already exist
files=".fna.gz"

# figure out which folders to exclude, e.g. unplaced_scaffolds, etc.

rsync -iPrLtm -f="+ *"$files"" -f="+ */" -f="- *" -f="- all_assembly_versions" ftp.ncbi.nlm.nih.gov::genomes/genbank/"$location" "$directory"

# --exclude-from=FILE     read exclude patterns from FILE
# --include-from=FILE     read include patterns from FILE
# -c, --checksum              skip based on checksum, not mod-time & size
# -i, --itemize-changes       output a change-summary for all updates
# -P, --partial, --progress   show progress and put partially downloaded files in a folder
# -r, --recursive             recurse into directories
# -L, --copy-links            transform symlink into referent file/dir
# -t, --times                 preserve modification times
# -m, --prune-empty-dirs      prune empty directory chains from file-list

mkdir -p "$directory"_local/
sudo find "$directory" -type f -exec cp -t ""$directory"_local/" -- {} +

