import pandas as pd
import os

# wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt

summary = '/home/truthling/MGGen/renamed/assembly_summary.txt'
genome_files = '/home/truthling/MGGen/Acinetobacter_nosocomialis/genome_files/'
df = pd.read_csv(summary, delimiter='\t', index_col=0)
copiedfiles = '/home/truthling/MGGen/renamed/'

for f in os.listdir(copiedfiles):
    if f.startswith('GCA'):
        id = (f.split('_')[0:2])
        id = ('_'.join(id))

    for row in df.index:
        if id == row:
            org_name = df.get_value(id, 'organism_name')
            strain = df.get_value(id, 'infraspecific_name').strip('strain=')
            assembly_level  = df.get_value(id, 'assembly_level')
            newname = '{}_{}_{}_{}.fasta'.format(id, org_name, strain, assembly_level)
            old = copiedfiles+f
            new = copiedfiles+newname
            os.rename(old, new)
