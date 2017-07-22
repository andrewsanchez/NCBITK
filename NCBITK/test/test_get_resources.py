from NCBITK import config
from NCBITK import get_resources

import unittest
import os
import re
import glob
import tempfile
import shutil
import pandas as pd


class TestGetResources(unittest.TestCase):
    def setUp(self):
        self.genbank_mirror = tempfile.mkdtemp(prefix='Genbank_')
        self.path_vars = config.instantiate_path_vars(self.genbank_mirror)
        self.info_dir, self.slurm, self.out, self.logger = self.path_vars
        self.assembly_summary = get_resources.get_assembly_summary(
            self.genbank_mirror, True)
        self.assertIsInstance(self.assembly_summary, pd.DataFrame)
        self.local_assembly_summary = pd.read_csv(
            'NCBITK/test/resources/assembly_summary.txt',
            sep="\t",
            index_col=0)
        self.path_assembly_summary = os.path.join(self.info_dir,
                                                  'assembly_summary.txt')

    def test_get_scientific_names(self):
        names = get_resources.get_scientific_names(self.genbank_mirror,
                                                   self.assembly_summary)
        self.assertIsInstance(names, pd.DataFrame)

    def test_update_assembly_summary(self):
        names = get_resources.get_scientific_names(self.genbank_mirror,
                                                   self.assembly_summary)
        updated_assembly_summary = get_resources.update_assembly_summary(
            self.assembly_summary, names)
        self.assertIsInstance(updated_assembly_summary, pd.DataFrame)

    def test_clean_up_assembly_summary(self):
        get_resources.clean_up_assembly_summary(
            self.local_assembly_summary)
        self.assertIsInstance(self.local_assembly_summary, pd.DataFrame)
        for i in self.local_assembly_summary.infraspecific_name:
            try:
                self.assertFalse(re.match('.*=.*', i) is True)
            except TypeError:
                continue

    def tearDown(self):
        shutil.rmtree(self.genbank_mirror)


if __name__ == '__main__':
    unittest.main()
