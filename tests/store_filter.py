import unittest
import pysam

from marker_alignments.store import AlignmentStore

class StoreFilter(unittest.TestCase):

    def assertStoreContent(self, alignment_store, expected):
        content = [t for t in alignment_store.query('select * from alignment')]
        self.assertEqual(content, expected)

    def test_filter_reads(self):
        r111 = ('taxon_1', 'marker_1','query_1', 1.0, 1.0)
        r112 = ('taxon_1', 'marker_1','query_2', 1.0, 1.0)
        r223 = ('taxon_2', 'marker_2','query_3', 1.0, 1.0)

        alignment_store = AlignmentStore()
        for r in [r111, r112, r223]:
          alignment_store.add_alignment(*r)
        self.assertStoreContent(alignment_store, [r111, r112, r223])

        alignment_store.modify_table_filter_taxa_on_num_markers_and_reads(min_num_markers = 1, min_num_reads = 2)
        self.assertStoreContent(alignment_store, [r111, r112])

    def test_filter_markers(self):
        r111 = ('taxon_1', 'marker_1','query_1', 1.0, 1.0)
        r122 = ('taxon_1', 'marker_2','query_2', 1.0, 1.0)
        r223 = ('taxon_2', 'marker_2','query_3', 1.0, 1.0)

        alignment_store = AlignmentStore()
        for r in [r111, r122, r223]:
          alignment_store.add_alignment(*r)

        alignment_store.modify_table_filter_taxa_on_num_markers_and_reads(min_num_markers = 2, min_num_reads = 0)
        self.assertStoreContent(alignment_store, [r111, r122])

    def test_filter_markers_best_matches_only(self):
        r111 = ('taxon_1', 'marker_1','query_1', 1.0, 1.0)
        r122 = ('taxon_1', 'marker_2','query_2', 0.5, 1.0)
        r222 = ('taxon_2', 'marker_2','query_2', 0.6, 1.0)

        alignment_store = AlignmentStore()
        for r in [r111, r122, r222]:
          alignment_store.add_alignment(*r)
        alignment_store.modify_table_filter_taxa_on_num_markers_and_reads(min_num_markers = 2, min_num_reads = 0)
        self.assertStoreContent(alignment_store, [])

    def test_filter_multiple_matches(self):
        r111 = ('taxon_1', 'marker_1','query_1', 1.0, 1.0)
        r122 = ('taxon_1', 'marker_2','query_2', 0.5, 1.0)
        r222 = ('taxon_2', 'marker_2','query_2', 0.6, 1.0)
        r333 = ('taxon_3', 'marker_3','query_3', 1.0, 1.0)

        alignment_store = AlignmentStore()
        for r in [r111, r122, r222, r333]:
          alignment_store.add_alignment(*r)

        # ask for all primary or unique matches
        alignment_store.modify_table_filter_taxa_on_multiple_matches(min_fraction_primary_matches = 1.0)
        # t1 has a secondary match and so gets filtered out
        self.assertStoreContent(alignment_store, [r222, r333])

if __name__ == '__main__':
    unittest.main()
