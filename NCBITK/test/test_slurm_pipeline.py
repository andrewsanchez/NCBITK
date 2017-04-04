import unittest
import os
import tempfile
import subprocess
import shutil
import pandas as pd

from NCBITK import curate


class TestSlurmPipeline(unittest.TestCase):

    def setUp(self):

        self.genbank_mirror = tempfile.mkdtemp()
        os.mkdir(os.path.join(self.genbank_mirror, ".info"))
        self.specification = 'NCBITK/test/specification.json'
        self.assembly_summary = pd.read_csv('NCBITK/test/resources/assembly_summary.txt', sep="\t", index_col=0)
        self.species_list = 'Acinetobacter_nosocomialis'

    def test_slurm_pipeline(self):

        self.cmd = "slurm-pipeline.py -s {} {}".format(self.specification, self.genbank_mirror)
        subprocess.Popen(self.cmd, shell=True).wait()
        self.assertTrue(os.listdir(self.genbank_mirror))

    def test_assess_genbank_mirror(self):

        genbank_assessment = curate.assess_genbank_mirror(self.genbank_mirror, self.assembly_summary, self.species_list)
        local_genomes, new_genomes, old_genomes, sketch_files, missing_sketch_files = genbank_assessment

        self.assertTrue(len(new_genomes) > len(local_genomes))
        self.assertTrue(len(local_genomes) == len(sketch_files))
        self.assertTrue(len(missing_sketch_files) > len(sketch_files))

    def tearDown(self):
        shutil.rmtree(self.genbank_mirror)
