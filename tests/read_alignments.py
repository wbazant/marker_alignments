import unittest
import pysam

import os 
dir_path = os.path.dirname(os.path.realpath(__file__))

import re
from marker_alignments.main import read_alignments

class ReadAlignments(unittest.TestCase):

    def test_example(self):
        sam = pysam.AlignmentFile(dir_path + "/data/example.sam")
        pattern_taxon = re.compile("^([^:]+):[^:]+$")
        pattern_marker = re.compile("^[^:]+:([^:]+)$")
        marker_to_taxon_id = { "taxon_1:marker_1" : "id_1"}
        alignment_store = read_alignments(sam, None, pattern_taxon, pattern_marker, marker_to_taxon_id)
        content = [t for t in alignment_store.query('select * from alignment')]
        self.assertEqual(content, [('taxon_1', 'marker_1', 'query_id', 1.0, 0.141026), ('taxon_2', 'marker_2', 'second_query_id', 0.973511, 0.02327), ('taxon_2', 'marker_3', 'third_query_id', 0.973511, 0.074),  ('taxon_2', 'marker_3', 'fourth_query_id', 0.973511, 0.074)])


if __name__ == '__main__':
    unittest.main()
