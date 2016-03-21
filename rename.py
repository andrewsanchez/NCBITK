import pandas as pd
import os

# wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt

summary = '/home/truthling/MGGen/assembly_summary.txt'
df = pd.read_csv(summary, delimiter='\t', index_col=0)
copiedfiles = '/home/truthling/MGGen/Acinetobacter_nosocomialis_local/'

df.update(df['infraspecific_name'][(df['infraspecific_name'].isnull()) & (df['isolate'].isnull())].fillna('NA'))
df.update(df['infraspecific_name'][(df['infraspecific_name'].isnull()) & (df['isolate'].notnull())].fillna(df['isolate']))
negate_list = ' '
df.organism_name.replace({'[^\w'+negate_list+']': ''}, regex=True, inplace=True)
df.infraspecific_name.replace({'[ =]': '_'}, regex=True, inplace=True)
df.assembly_level.replace({' ': '_'}, regex=True, inplace=True)

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
