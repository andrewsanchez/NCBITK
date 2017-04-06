#!/bin/bash

species='Acinetobacter_nosocomialis'
python ~/projects/NCBITK/run.py --update "$@" --use_local --species $species
# python ~/projects/NCBITK/NCBITK/filtering/run.py -p "$@" --use_local --species $species
