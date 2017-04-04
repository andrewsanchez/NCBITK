#!/bin/bash

pwd
mkdir -p "$@"/.info
cp NCBITK/test/resources/assembly_summary.txt "$@"/.info/assembly_summary.txt
python ~/projects/NCBITK/run.py "$@" --use_local --species Acinetobacter_nosocomialis
ls -alR "$@"
