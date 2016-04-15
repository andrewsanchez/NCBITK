import os

emptyDirs = []

for root, dirs, files in os.walk('genbank_bacteria'):
    for name in dirs:
        name = os.path.join(root, name)
        if not os.listdir(name):
            print(name + ' is empty.')
            emptyDirs.append(name)
print('Empty directories:  ' + str(len(emptyDirs)))
