import os, sys
import pandas as pd

renametarget = sys.argv[1]

# clean up assembly_summary.txt
df = pd.read_csv('assembly_summary.txt', delimiter='\t', index_col=0)
df.update(df['infraspecific_name'][(df['infraspecific_name'].isnull()) & (df['isolate'].isnull())].fillna('NA'))
df.update(df['infraspecific_name'][(df['infraspecific_name'].isnull()) & (df['isolate'].notnull())].fillna(df['isolate']))
df.assembly_level.replace({' ': '_'}, regex=True, inplace=True)
df.organism_name.replace({' = ': '_'}, regex=True, inplace=True)
df.organism_name.replace({' ': '_'}, regex=True, inplace=True)
df.organism_name.replace({'[\W]': '_'}, regex=True, inplace=True)
df.infraspecific_name.replace({' = ': '_'}, regex=True, inplace=True)
df.infraspecific_name.replace({'[\W]': '_'}, regex=True, inplace=True)

# remove duplicate strings during renaming
def rmDuplicates(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]

for root, dirs, files in os.walk(renametarget):
    for f in files:
        id = (f.split('_')[0:2])
        id = ('_'.join(id))
        if id in df.index:
           #f = os.path.join(root, f)
           #print(f)
            org_name = df.get_value(id, 'organism_name')
            strain = df.get_value(id, 'infraspecific_name').strip('strain=')
            assembly_level  = df.get_value(id, 'assembly_level')
            newname = '{}_{}{}_{}.fna'.format(id, org_name, strain, assembly_level)
            newname = newname.split('_')
            newname = rmDuplicates(newname)
            newname = '_'.join(newname)
            newname = newname.replace('subsp_', '')
            newname = newname.replace('str_', '')
            newname = newname.replace('strain_', '')
            old = os.path.join(root, f)
            new = os.path.join(root, newname)
            os.rename(old, new)
