#!/usr/bin/env python

import sys
import subprocess
import os

organism = sys.argv[1]
local_mirror = sys.argv[2]
ftp_directory = 'bacteria/' + organism + '/latest_assembly_versions/'
just_fastas = local_mirror + '_local'

subprocess.run(['rsync',
                '-iPrLtmn',
                '-f="- **unplaced_scaffolds**"',
                '-f="+ *.fna.gz"',
                '-f="+ */"',
                '-f="- *"',
                'ftp.ncbi.nlm.nih.gov::genomes/genbank/' + ftp_directory,
                local_mirror,])

if os.path.isdir(just_fastas):
    subprocess.run(['find', organism, '-type', 'f',
                '-exec', 'cp',
                '-t', just_fastas,
                '-- {}', '+'])

else:
    os.mkdir(just_fastas)
    subprocess.run(['find', organism, '-type', 'f',
                '-exec', 'cp',
                '-t', just_fastas,
                '-- {}', '+'])

# sudo find "$ftpmirror" -type f -exec cp -t "$justfastas" -- {} +
