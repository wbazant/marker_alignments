import unittest

import sys
import os 
dir_path = os.path.dirname(os.path.realpath(__file__))
input_path = dir_path + "/data/marker_to_taxon.tsv"

from marker_alignments.main import read_marker_to_taxon

class MarkerToTaxon(unittest.TestCase):
    def test_read_marker_to_taxon(self):
        (marker_to_taxon, taxon_to_num_markers) = read_marker_to_taxon(input_path)
        self.assertEqual(marker_to_taxon , {"m1": "t1", "m21": "t2", "m22": "t2"}, "Marker to taxon")
        self.assertEqual(taxon_to_num_markers,{"t1": 1, "t2": 2} , "Taxon to num markers")

if __name__ == '__main__':
    unittest.main()

