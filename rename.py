import os, re, argparse
import pandas as pd

# remove duplicate strings during renaming
def rmduplicates(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]

def rename(renametarget):
    # clean up assembly_summary.txt
    df = pd.read_csv('assembly_summary.txt', delimiter='\t', index_col=0)
    df.update(df['infraspecific_name'][(df['infraspecific_name'].isnull()) & (df['isolate'].isnull())].fillna('NA'))
    df.update(df['infraspecific_name'][(df['infraspecific_name'].isnull()) & (df['isolate'].notnull())].fillna(df['isolate']))
    df.assembly_level.replace({' ': '_'}, regex=True, inplace=True)
    df.organism_name.replace({' ': '_'}, regex=True, inplace=True)
    df.organism_name.replace({'[\W]': '_'}, regex=True, inplace=True)
    df.infraspecific_name.replace({'[\W]': '_'}, regex=True, inplace=True)

    for root, dirs, files in os.walk(renametarget):
        for f in files:
            id = (f.split('_')[0:2])
            id = ('_'.join(id))
            if id in df.index:
                org_name = df.get_value(id, 'organism_name')
                strain = df.get_value(id, 'infraspecific_name')
                assembly_level  = df.get_value(id, 'assembly_level')
                newname = '{}_{}_{}_{}.fna'.format(id, org_name, strain, assembly_level)
                rmwords = re.compile( r'((?<=_)(sp|sub|substr|subsp|str|strain)(?=_))' )
                newname = rmwords.sub('_', newname)
                newname = re.sub(r'_+', '_', newname )
                newname = newname.split('_')
                newname = rmduplicates(newname)
                newname = '_'.join(newname)
                print(newname)
                old = os.path.join(root, f)
                new = os.path.join(root, newname)
                os.rename(old, new)

def Main():
    parser = argparse.ArgumentParser()
    parser.add_argument('rename_target', help = 'The folder whose contents will be renamed', type=str)
    args = parser.parse_args()

    rename(args.rename_target)

if __name__ == '__main__':
    Main()
