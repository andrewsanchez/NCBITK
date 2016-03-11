#!/bin/bash

# Run as sudo
# include="rsync_include.txt"
# exclude="rsync_exclude.txt"
organism="Acinetobacter_nosocomialis" # set to "*" for all bacteria
location="bacteria/"$organism"/latest_assembly_versions/"
directory="/home/truthling/MGGen/rsync_Acinetobacter_nosocomialis/" # will be created if doesn't already exist
files=".fna.gz"
renamescript="/home/truthling/MGGen/NCBI_tools/rename.python"

# figure out which folders to exclude, e.g. unplaced_scaffolds, etc.

rsync -iPrLtm -f="+ *"$files"" -f="+ */" -f="- *" ftp.ncbi.nlm.nih.gov::genomes/genbank/"$location" "$directory"

# these two below produce the same results as the one above
# rsync -iPrLtm -f="+ *"$files"" -f="+ */" -f="- *" ftp.ncbi.nlm.nih.gov::genomes/genbank/"$location" $directory 
# rsync -iPrLtm -f="+ *"$files"" -f="+ */" -f="- *" -f="- all_assembly_versions" ftp.ncbi.nlm.nih.gov::genomes/genbank/"$location" "$directory"

# cp only the files to separate directory for renaming
# Make sure this is only running once, after all files are downloaded.

find "$directory" -type f -exec cp -t '/home/truthling/MGGen/renamed/' -- {} +

python3 "$renamescript"

# add line for decompressing files

# --exclude-from=FILE     read exclude patterns from FILE
# --include-from=FILE     read include patterns from FILE
# -c, --checksum              skip based on checksum, not mod-time & size
# -i, --itemize-changes       output a change-summary for all updates
# -P, --partial, --progress   show progress and put partially downloaded files in a folder
# -r, --recursive             recurse into directories
# -L, --copy-links            transform symlink into referent file/dir
# -t, --times                 preserve modification times
# -m, --prune-empty-dirs      prune empty directory chains from file-list
