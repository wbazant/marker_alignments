import unittest

import pysam

import os
dir_path = os.path.dirname(os.path.realpath(__file__))

import re
from marker_alignments.main import get_output, output_type_options as output_types
from marker_alignments.store import AlignmentStore

num_reads = 1000000
def output_types_with_column(value_column):
    result = []
    for output_type in output_types:
        header, lines = get_output(AlignmentStore(), output_type, num_reads)
        if value_column in header:
            result.append(output_type)
    return result


def get_as_dictionary(alignments, output_type, value_column):
    alignment_store = AlignmentStore()
    for t in alignments:
        alignment_store.add_alignment(*t)
    header, lines = get_output(alignment_store, output_type, num_reads)

    tix = header.index("taxon")
    mix = header.index("marker") if "marker" in header else None
    vix = header.index(value_column)

    result = {}
    for l in lines:
        if "marker" in header:
            if l[tix] not in result:
                result[l[tix]] = {}
            if l[mix] not in result[l[tix]]:
                result[l[tix]][l[mix]] = {}
            result[l[tix]][l[mix]] = l[vix]
        else:
            if l[tix] not in result:
                result[l[tix]] = {}
            result[l[tix]] = l[vix]
    return result

def test_column(self, alignments, value_column, expected):
    for output_type in output_types_with_column(value_column):
        with self.subTest(output_type):
            self.assertEqual(get_as_dictionary(alignments, output_type, value_column), expected)

# rxyz = (taxon_x, marker_y, query_z, identity, coverage)
class Store(unittest.TestCase):
    def test_counting(self):
        r111 = ('taxon_1', 'marker_1','query_1', 1.0, 1.0)

        test_column(self, [r111], "marker_read_count", {"taxon_1" : {"marker_1" : 1}})
        test_column(self, [r111], "marker_cpm", {"taxon_1" : {"marker_1" : 1.0}})
        test_column(self, [r111], "cpm", {"taxon_1" : 1.0})
        test_column(self, [r111], "taxon_num_reads", {"taxon_1" : 1.0})
        test_column(self, [r111], "taxon_num_markers", {"taxon_1" : 1.0})

        r112 = ('taxon_1', 'marker_1','query_2', 1.0, 1.0)
        test_column(self, [r111, r112], "marker_read_count", {"taxon_1" : {"marker_1" : 2}})
        test_column(self, [r111, r112], "marker_cpm", {"taxon_1" : {"marker_1" : 2}})
        test_column(self, [r111, r112], "cpm", {"taxon_1" : 2.0})
        test_column(self, [r111, r112], "taxon_num_reads", {"taxon_1" : 2.0})
        test_column(self, [r111, r112], "taxon_num_markers", {"taxon_1" : 1.0})

        r223 = ('taxon_2', 'marker_2','query_3', 1.0, 1.0)
        test_column(self, [r111, r112, r223], "marker_read_count", {"taxon_1" : {"marker_1" : 2}, "taxon_2": {"marker_2": 1}})
        test_column(self, [r111, r112, r223], "marker_cpm", {"taxon_1" : {"marker_1" : 2}, "taxon_2": {"marker_2": 1}})
        test_column(self, [r111, r112, r223], "cpm", {"taxon_1" : 2.0, "taxon_2": 1.0})
        test_column(self, [r111, r112, r223], "taxon_num_reads", {"taxon_1" : 2.0, "taxon_2": 1.0})
        test_column(self, [r111, r112, r223], "taxon_num_markers", {"taxon_1" : 1.0, "taxon_2": 1.0})

        # match split over two markers
        r114 = ('taxon_1', 'marker_1','query_4', 1.0, 1.0)
        r224 = ('taxon_2', 'marker_2','query_4', 1.0, 1.0)
        test_column(self, [r111, r112, r223, r114, r224], "marker_read_count", {"taxon_1" : {"marker_1" : 2.5}, "taxon_2": {"marker_2": 1.5}})
        test_column(self, [r111, r112, r223, r114, r224], "marker_cpm", {"taxon_1" : {"marker_1" : 2.5}, "taxon_2": {"marker_2": 1.5}})
        test_column(self, [r111, r112, r223, r114, r224], "cpm", {"taxon_1" : 2.5, "taxon_2": 1.5})
        test_column(self, [r111, r112, r223, r114, r224], "taxon_num_reads", {"taxon_1" : 2.5, "taxon_2": 1.5})
        test_column(self, [r111, r112, r223, r114, r224], "taxon_num_markers", {"taxon_1" : 1, "taxon_2": 1})

        # additional match to a different marker
        r155 = ('taxon_1', 'marker_3','query_5', 1.0, 1.0)
        test_column(self, [r111, r155], "marker_read_count", {"taxon_1" : {"marker_1" : 1, "marker_3": 1}})
        test_column(self, [r111, r155], "marker_cpm", {"taxon_1" : {"marker_1" : 1, "marker_3": 1}})
        test_column(self, [r111, r155], "cpm", {"taxon_1" : 1.0})
        test_column(self, [r111, r155], "taxon_num_reads", {"taxon_1" : 2})
        test_column(self, [r111, r155], "taxon_num_markers", {"taxon_1" : 2})

        # an uneven match split over two markers
        r116 = ('taxon_1', 'marker_1','query_6', 0.9, 1.0)
        r226 = ('taxon_2', 'marker_2','query_6', 0.7, 1.0)
        test_column(self, [r111, r112, r223, r116, r226], "marker_read_count", {"taxon_1" : {"marker_1" : 2.623076923076923}, "taxon_2": {"marker_2": 1.376923076923077}})
        test_column(self, [r111, r112, r223, r116, r226], "marker_cpm", {"taxon_1" : {"marker_1" : 2.623076923076923}, "taxon_2": {"marker_2": 1.376923076923077}})
        test_column(self, [r111, r112, r223, r116, r226], "cpm", {"taxon_1" : 2.623076923076923, "taxon_2": 1.376923076923077})

        # a shorter read that aligns only to 0.25 of a marker
        r117 = ('taxon_1', 'marker_1','query_7', 1.0, 0.25)
        test_column(self, [r117], "marker_read_count", {"taxon_1" : {"marker_1" : 1}})
        test_column(self, [r117], "marker_cpm", {"taxon_1" : {"marker_1" : 0.25}})
        test_column(self, [r117], "cpm", {"taxon_1" : 0.25})
        test_column(self, [r117], "taxon_num_reads", {"taxon_1" : 1.0})
        test_column(self, [r117], "taxon_num_markers", {"taxon_1" : 1.0})

    def test_identity(self):
        r111 = ('taxon_1', 'marker_1','query_1', 0.5, 1.0)
        r112 = ('taxon_1', 'marker_1','query_2', 1.0, 1.0)
        test_column(self, [r111], "marker_avg_identity", {"taxon_1" : {"marker_1" : 0.5}})
        test_column(self, [r111, r112], "marker_avg_identity", {"taxon_1" : {"marker_1" : 0.75}})

if __name__ == '__main__':
    unittest.main()

