import sys
import subprocess

#organism="$1"
#ftplocation="bacteria/"$organism"/latest_assembly_versions/"
#ftpmirror="./"$organism""
#justfastas="$ftpmirror"_local
#mkdir -p "$justfastas"

organism = sys.argv[1]
local_mirrior = sys.argv[2]
ftp_directory = 'bacteria/' + organism + '/latest_assembly_versions/'

if organism = '*':
subprocess.run(['rsync',
                '-iPrLtm',
                '-f="- **unplaced_scaffolds**"',
                '-f="+ *.fna.gz"',
                '-f="+ */"',
                '-f="- *"',
                'ftp.ncbi.nlm.nih.gov::genomes/genbank/' + ftp_directory,
                '"$ftpmirror"',
                shell=True
])

rsync -iPrLtm -f="- **unplaced_scaffolds**" -f="+ *.fna.gz" -f="+ */" -f="- *" ftp.ncbi.nlm.nih.gov::genomes/genbank/"$ftplocation" "$ftpmirror"

sudo find "$ftpmirror" -type f -exec cp -t "$justfastas" -- {} +
