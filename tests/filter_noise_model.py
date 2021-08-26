import unittest

import logging
from marker_alignments.filter.noise_model import fit_noise_model

def t(taxon_to_num_markers):
    return fit_noise_model(taxon_to_num_markers, 10000)

class NoiseModel(unittest.TestCase):
    def test_zeros_at_the_end(self):
        self.assertEqual(t({0:100,1:0}), [(0,1),(1,1)])

    def test_ones(self):
        r = t({0:95, 1: 5})
        (n,p) = r[1]
        self.assertGreater(p, 0.05)

    def test_zeros_at_the_end_dont_change_result(self):
        r1 = t({0:95, 1: 5})
        (n11,p11) = r1[1]
        r2 = t({0:95, 1: 5, 2:0})
        (n21,p21) = r2[1]
        self.assertEqual(p11,p21)

    def test_large_numbers_at_the_end_make_earlier_numbers_more_probable(self):
        r1 = t({0:95, 1: 5})
        (n11,p11) = r1[1]
        r2 = t({0:95, 1: 5, 2:10})
        (n21,p21) = r2[1]
        self.assertGreater(p21, p11)

    def test_large_number_at_the_end_is_improbable(self):
        r = t({0:75, 1: 5, 2:10})
        (n,p) = r[2]
        self.assertLess(p, 0.05)

    def test_small_number_at_the_end_is_probable(self):
        r = t({0:944, 1: 50, 2:3})
        (n,p) = r[2]
        self.assertGreater(p, 0.05)

    def test_medium_number_at_the_end_is_vaguely_improbable(self):
        r = t({0:941, 1: 50, 2:3, 3:1})
        (n,p) = r[3]
        self.assertLess(p, 0.05)
        self.assertGreater(p, 0.01)

    def skip_test_typical_dataset(self):
        r = t({0: 2250, 1:1117,2:466,3:137,4:35,5:7,6:3,7:1,25:1})
        self.assertEqual(r,[])

if __name__ == '__main__':
    unittest.main()
