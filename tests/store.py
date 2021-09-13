import unittest
import pysam

import os 
dir_path = os.path.dirname(os.path.realpath(__file__))

import re
from marker_alignments.main import read_alignments
from marker_alignments.store import AlignmentStore

num_reads = 100000
def outputs(alignment_store):
    return {
      "marker_all": sorted([t for t in alignment_store.as_marker_all(num_reads)]),
      "marker_coverage": sorted([t for t in alignment_store.as_marker_coverage()]),
      "marker_cpm": sorted([t for t in alignment_store.as_marker_cpm(num_reads)]),
      "marker_read_count": sorted([t for t in alignment_store.as_marker_read_count()]),
      "taxon_all": sorted([t for t in alignment_store.as_taxon_all(num_reads)]),
      "taxon_coverage": sorted([t for t in alignment_store.as_taxon_coverage()]),
      "taxon_cpm": sorted([t for t in alignment_store.as_taxon_cpm(num_reads)]),
      "taxon_read_and_marker_count": sorted([t for t in alignment_store.as_taxon_read_and_marker_count()]),
    }

class MarkerScores(unittest.TestCase):

    def results_as_expected(self, alignments, expected_scores={}):
        self.maxDiff=None
        alignment_store = AlignmentStore()
        for t in alignments:
            alignment_store.add_alignment(*t)

        self.assertEqual(outputs(alignment_store),expected_scores) 

    def test_one_read(self):
        self.results_as_expected(
            alignments = [
                ('taxon_1', 'marker_1','query_1', 1.0, 0.111),
            ],
            expected_scores = {'marker_all': [('taxon_1', 'marker_1', 0.111, 1.1099999999999999, 1.0)],
              'marker_coverage': [('taxon_1', 'marker_1', 0.111)],
              'marker_cpm': [('taxon_1', 'marker_1', 1.1099999999999999)],
              'marker_read_count': [('taxon_1', 'marker_1', 1.0)],
              'taxon_all': [('taxon_1', 0.111, 1.1099999999999999, 1.0, 1, 1.0)],
              'taxon_coverage': [('taxon_1', 0.111)],
              'taxon_cpm': [('taxon_1', 1.1099999999999999)],
              'taxon_read_and_marker_count': [('taxon_1', 1.0, 1, 1.0)]}
        )

    def test_one_read_two_markers(self):
        self.results_as_expected(
            alignments = [
                ('taxon_1', 'marker_1','query_1', 0.11, 0.444),
                ('taxon_1', 'marker_2','query_1', 0.33, 0.444),
            ],
            expected_scores = {'marker_all': [('taxon_1', 'marker_1', 0.111, 1.1099999999999999, 0.25),
                             ('taxon_1', 'marker_2', 0.333, 3.3300000000000005, 0.75)],
              'marker_coverage': [('taxon_1', 'marker_1', 0.111),
                                  ('taxon_1', 'marker_2', 0.333)],
              'marker_cpm': [('taxon_1', 'marker_1', 1.1099999999999999),
                             ('taxon_1', 'marker_2', 3.3300000000000005)],
              'marker_read_count': [('taxon_1', 'marker_1', 0.25),
                                    ('taxon_1', 'marker_2', 0.75)],
              'taxon_all': [('taxon_1', 0.222, 2.2199999999999998, 1.0, 2, 0.75)],
              'taxon_coverage': [('taxon_1', 0.222)],
              'taxon_cpm': [('taxon_1', 2.2199999999999998)],
              'taxon_read_and_marker_count': [('taxon_1', 1.0, 2, 0.75)]}
        )

    def test_one_read_two_markers_with_extra_read(self):
        self.results_as_expected(
            alignments = [
                ('taxon_1', 'marker_1','query_1', 0.11, 0.444),
                ('taxon_1', 'marker_2','query_1', 0.33, 0.444),
                ('taxon_1', 'marker_1','query_2', 1, 0.5),
            ],
            expected_scores = {'marker_all': [('taxon_1', 'marker_1', 0.611, 6.11, 1.25),
                             ('taxon_1', 'marker_2', 0.333, 3.3300000000000005, 0.75)],
              'marker_coverage': [('taxon_1', 'marker_1', 0.611),
                                  ('taxon_1', 'marker_2', 0.333)],
              'marker_cpm': [('taxon_1', 'marker_1', 6.11),
                             ('taxon_1', 'marker_2', 3.3300000000000005)],
              'marker_read_count': [('taxon_1', 'marker_1', 1.25),
                                    ('taxon_1', 'marker_2', 0.75)],
              'taxon_all': [('taxon_1', 0.472, 4.72, 2.0, 2, 1.25)],
              'taxon_coverage': [('taxon_1', 0.472)],
              'taxon_cpm': [('taxon_1', 4.72)],
              'taxon_read_and_marker_count': [('taxon_1', 2.0, 2, 1.25)]}
        )

    def test_one_read_weight_does_not_matter(self):
        self.results_as_expected(
            alignments = [
                ('taxon_1', 'marker_1','query_1', 0.42, 0.111),
            ],
            expected_scores = {'marker_all': [('taxon_1', 'marker_1', 0.111, 1.1099999999999999, 1.0)],
              'marker_coverage': [('taxon_1', 'marker_1', 0.111)],
              'marker_cpm': [('taxon_1', 'marker_1', 1.1099999999999999)],
              'marker_read_count': [('taxon_1', 'marker_1', 1.0)],
              'taxon_all': [('taxon_1', 0.111, 1.1099999999999999, 1.0, 1, 1.0)],
              'taxon_coverage': [('taxon_1', 0.111)],
              'taxon_cpm': [('taxon_1', 1.1099999999999999)],
              'taxon_read_and_marker_count': [('taxon_1', 1.0, 1, 1.0)]}
        )

    def test_two_reads_one_marker(self):
        self.results_as_expected(
            alignments = [
                ('taxon_1', 'marker_1','query_1', 1.0, 0.111),
                ('taxon_1', 'marker_1','query_2', 1.0, 0.112),
            ],
            expected_scores = {'marker_all': [('taxon_1', 'marker_1', 0.223, 2.2300000000000004, 2.0)],
              'marker_coverage': [('taxon_1', 'marker_1', 0.223)],
              'marker_cpm': [('taxon_1', 'marker_1', 2.2300000000000004)],
              'marker_read_count': [('taxon_1', 'marker_1', 2.0)],
              'taxon_all': [('taxon_1', 0.223, 2.2300000000000004, 2.0, 1, 2.0)],
              'taxon_coverage': [('taxon_1', 0.223)],
              'taxon_cpm': [('taxon_1', 2.2300000000000004)],
              'taxon_read_and_marker_count': [('taxon_1', 2.0, 1, 2.0)]}
        )

    def test_two_reads_two_markers(self):
        self.results_as_expected(
            alignments = [
                ('taxon_1', 'marker_1','query_1', 1.0, 0.111),
                ('taxon_1', 'marker_2','query_2', 1.0, 0.122),
            ],
            expected_scores = {'marker_all': [('taxon_1', 'marker_1', 0.111, 1.1099999999999999, 1.0),
                             ('taxon_1', 'marker_2', 0.122, 1.22, 1.0)],
              'marker_coverage': [('taxon_1', 'marker_1', 0.111),
                                  ('taxon_1', 'marker_2', 0.122)],
              'marker_cpm': [('taxon_1', 'marker_1', 1.1099999999999999),
                             ('taxon_1', 'marker_2', 1.22)],
              'marker_read_count': [('taxon_1', 'marker_1', 1.0),
                                    ('taxon_1', 'marker_2', 1.0)],
              'taxon_all': [('taxon_1', 0.11649999999999999, 1.1649999999999998, 2.0, 2, 1.0)],
              'taxon_coverage': [('taxon_1', 0.11649999999999999)],
              'taxon_cpm': [('taxon_1', 1.1649999999999998)],
              'taxon_read_and_marker_count': [('taxon_1', 2.0, 2, 1.0)]}
        )

    def test_two_reads_two_taxons(self):
        self.results_as_expected(
            alignments = [
                ('taxon_1', 'marker_1','query_1', 1.0, 0.111),
                ('taxon_2', 'marker_2','query_2', 1.0, 0.222),
            ],
            expected_scores = {'marker_all': [('taxon_1', 'marker_1', 0.111, 1.1099999999999999, 1.0),
                             ('taxon_2', 'marker_2', 0.222, 2.2199999999999998, 1.0)],
              'marker_coverage': [('taxon_1', 'marker_1', 0.111),
                                  ('taxon_2', 'marker_2', 0.222)],
              'marker_cpm': [('taxon_1', 'marker_1', 1.1099999999999999),
                             ('taxon_2', 'marker_2', 2.2199999999999998)],
              'marker_read_count': [('taxon_1', 'marker_1', 1.0),
                                    ('taxon_2', 'marker_2', 1.0)],
              'taxon_all': [('taxon_1', 0.111, 1.1099999999999999, 1.0, 1, 1.0),
                            ('taxon_2', 0.222, 2.2199999999999998, 1.0, 1, 1.0)],
              'taxon_coverage': [('taxon_1', 0.111), ('taxon_2', 0.222)],
              'taxon_cpm': [('taxon_1', 1.1099999999999999),
                            ('taxon_2', 2.2199999999999998)],
              'taxon_read_and_marker_count': [('taxon_1', 1.0, 1, 1.0), ('taxon_2', 1.0, 1, 1.0)]}
        )


if __name__ == '__main__':
    unittest.main()

