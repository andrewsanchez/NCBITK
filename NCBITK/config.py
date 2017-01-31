#!/usr/bin/env python

import os

def instantiate_path_vars(genbank_mirror):

    info_dir = os.path.join(genbank_mirror, ".info")
    slurm = os.path.join(info_dir, "slurm")
    out = os.path.join(slurm, "out")

    return info_dir, slurm, out
