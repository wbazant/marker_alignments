import unittest
import pysam

from marker_alignments.store import AlignmentStore

class StoreTransformTaxa(unittest.TestCase):

    def assertStoreContent(self, alignment_store, expected):
        content = [t[0] for t in alignment_store.query('select distinct taxon from alignment')]
        self.assertEqual(set(content), set(expected))

    def test_one(self):
        r111 = ('taxon_1', 'marker_1','query_1', 0.7, 1.0)
        
        clusters = [('taxon_1',)]
        alignment_store = AlignmentStore()
        for r in [r111]:
          alignment_store.add_alignment(*r)

        alignment_store._store_taxon_clusters(clusters)

        alignment_store.modify_table_transform_taxa_on_thresholds_and_clusters(
                threshold_identity = 0.4,
                min_num_taxa_below_identity = 0,
                min_num_markers_below_identity = 0,
                min_num_reads_below_identity =0)

        self.assertStoreContent(alignment_store, ['taxon_1'])

    def test_one_cluster(self):
        r111 = ('taxon_1', 'marker_1','query_1', 0.7, 1.0)
        r222 = ('taxon_2', 'marker_2','query_2', 0.7, 1.0)
        
        clusters = [('taxon_1','taxon_2')]
        alignment_store = AlignmentStore()
        for r in [r111, r222]:
          alignment_store.add_alignment(*r)

        alignment_store._store_taxon_clusters(clusters)

        alignment_store.modify_table_transform_taxa_on_thresholds_and_clusters(
                threshold_identity = 0.4,
                min_num_taxa_below_identity = 0,
                min_num_markers_below_identity = 0,
                min_num_reads_below_identity =0)

        self.assertStoreContent(alignment_store, ['taxon_1', 'taxon_2'])

    def test_one_cluster_filter_taxa_below_threshold(self):
        r111 = ('taxon_1', 'marker_1','query_1', 0.7, 1.0)
        r222 = ('taxon_2', 'marker_2','query_2', 0.3, 1.0)
        
        clusters = [('taxon_1','taxon_2')]
        alignment_store = AlignmentStore()
        for r in [r111, r222]:
          alignment_store.add_alignment(*r)

        alignment_store._store_taxon_clusters(clusters)

        alignment_store.modify_table_transform_taxa_on_thresholds_and_clusters(
                threshold_identity = 0.4,
                min_num_taxa_below_identity = 0,
                min_num_markers_below_identity = 0,
                min_num_reads_below_identity =0)

        self.assertStoreContent(alignment_store, ['taxon_1'])

    def test_one_cluster_merge_taxa_below_threshold(self):
        r111 = ('taxon_1', 'marker_1','query_1', 0.7, 1.0)
        r222 = ('taxon_2', 'marker_2','query_2', 0.7, 1.0)
        
        clusters = [('taxon_1','taxon_2')]
        alignment_store = AlignmentStore()
        for r in [r111, r222]:
          alignment_store.add_alignment(*r)

        alignment_store._store_taxon_clusters(clusters)

        alignment_store.modify_table_transform_taxa_on_thresholds_and_clusters(
                threshold_identity = 1,
                min_num_taxa_below_identity = 1,
                min_num_markers_below_identity = 0,
                min_num_reads_below_identity =0)

        self.assertStoreContent(alignment_store, ['?taxon_1,taxon_2'])

    def test_one_cluster_merge_taxa_below_threshold_2(self):
        r111 = ('taxon_1', 'marker_1','query_1', 1, 1.0)
        r222 = ('taxon_2', 'marker_2','query_2', 0.7, 1.0)
        r333 = ('taxon_3', 'marker_3','query_3', 0.7, 1.0)
        
        clusters = [('taxon_1', 'taxon_2', 'taxon_3')]
        alignment_store = AlignmentStore()
        for r in [r111, r222, r333]:
          alignment_store.add_alignment(*r)

        alignment_store._store_taxon_clusters(clusters)

        alignment_store.modify_table_transform_taxa_on_thresholds_and_clusters(
                threshold_identity = 1,
                min_num_taxa_below_identity = 1,
                min_num_markers_below_identity = 0,
                min_num_reads_below_identity =0)

        self.assertStoreContent(alignment_store, ['taxon_1'])

if __name__ == '__main__':
    unittest.main()
