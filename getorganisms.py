#!/usr/bin/env python3

def getorganisms(organismfile):
    dirs = [] 
    with open(organismfile, 'r') as getorganism:
        org = getorganism.readline()
        dirs.append(org)

    for org in dirs:
        print(org)

if __name__ == "__main__":
    import sys
    getorganisms(sys.argv[1])
