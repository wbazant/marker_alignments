import unittest
import pysam

import os 
dir_path = os.path.dirname(os.path.realpath(__file__))

import re
from marker_alignments.main import read_alignments

r1=('id_1|taxon_1', 'marker_1', 'query_id', 0.948052, 0.141026)
r2=('taxon_2', 'marker_2', 'second_query_id', 0.933333, 0.02327)
r3=('taxon_2', 'marker_3', 'third_query_id', 0.933333, 0.074)
r4=('taxon_2', 'marker_3', 'fourth_query_id', 0.933333, 0.074)
r5=('taxon_3', 'marker_1', 'fifth_query_id', 0.95, 0.02)

pattern_taxon = re.compile("^([^:]+):[^:]+$")
pattern_marker = re.compile("^[^:]+:([^:]+)$")
marker_to_taxon_id = { "taxon_1:marker_1" : "id_1"}

class ReadAlignments(unittest.TestCase):

    def assertStoreContent(self, alignment_store, expected):
        content = [t for t in alignment_store.query('select * from alignment')]
        self.assertEqual(content, expected)

    def test_example(self):
        sam = pysam.AlignmentFile(dir_path + "/data/example.sam")
        alignment_store = read_alignments(sam, None, pattern_taxon, pattern_marker, marker_to_taxon_id, 0,0,0)
        self.assertStoreContent(alignment_store, [r1,r2,r3,r4,r5])

    def test_filter_mapq(self):
        sam = pysam.AlignmentFile(dir_path + "/data/example.sam")
        alignment_store = read_alignments(sam, None, pattern_taxon, pattern_marker, marker_to_taxon_id, min_mapq = 10, min_query_length = 0, min_match_identity = 0)
        self.assertStoreContent(alignment_store, [r1,r2,r3,r4])

    def test_filter_query_length(self):
        sam = pysam.AlignmentFile(dir_path + "/data/example.sam")
        alignment_store = read_alignments(sam, None, pattern_taxon, pattern_marker, marker_to_taxon_id, min_mapq = 0, min_query_length = 25, min_match_identity = 0)
        self.assertStoreContent(alignment_store, [r1,r2,r3,r4])

    def test_filter_match_identity(self):
        sam = pysam.AlignmentFile(dir_path + "/data/example.sam")
        alignment_store = read_alignments(sam, None, pattern_taxon, pattern_marker, marker_to_taxon_id,  min_mapq = 0, min_query_length = 0, min_match_identity = 0.94 )
        self.assertStoreContent(alignment_store, [r1,r5])




if __name__ == '__main__':
    unittest.main()
