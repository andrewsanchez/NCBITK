#!/usr/bin/env python

import os
import time

def instantiate_path_vars(genbank_mirror):

    info_dir = os.path.join(genbank_mirror, ".info")
    slurm = os.path.join(info_dir, "slurm")
    out = os.path.join(slurm, "out")
    ymd = time.strftime("%y.%m.%d")
    log_file = os.path.join(info_dir, 'log_{}.txt'.format(ymd))

    return info_dir, slurm, out, log_file
