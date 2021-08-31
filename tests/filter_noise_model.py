import unittest

import logging
#logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('tests/filter_noise_model')
from marker_alignments.filter.noise_model import fit_noise_model, log_likelihood, cutoff_fit_for_noise_model


class LogLikelihood(unittest.TestCase):
    def test_decreasing_numbers_are_more_likely(self):
        ll = log_likelihood([1000, 500, 0], 10000)
        ll2 = log_likelihood([1000, 500, 100], 10000)
        ll3 = log_likelihood([1000, 500, 500], 10000)
        self.assertLess(ll, ll2)
        self.assertGreater(ll2, ll3)

class PickCutoff(unittest.TestCase):
    def test_null_case(self):
        cutoff = cutoff_fit_for_noise_model({0: 2250, 1:1117}, 6000, logger)
        self.assertEqual(cutoff, 2)

    def test_case_small(self):
        cutoff = cutoff_fit_for_noise_model({0: 4000, 1:21,2:1 }, 6000, logger)
        self.assertEqual(cutoff, 2)

    def test_case_prawn(self):
        cutoff = cutoff_fit_for_noise_model({0: 2250, 1:1117,2:466,3:137,4:35,5:7,6:3,7:1,25:1}, 6000, logger)
        self.assertEqual(cutoff, 8)

def noise_model(taxon_to_num_markers):
    return fit_noise_model(taxon_to_num_markers, 10000, logger)

class NoiseModel(unittest.TestCase):
    def test_zeros_at_the_end(self):
        self.assertEqual(noise_model({0:100,1:0}), [(0,1),(1,1)])

    def test_ones(self):
        r = noise_model({0:95, 1: 5})
        (n,p) = r[1]
        self.assertGreater(p, 0.05)

    def test_zeros_at_the_end_dont_change_result(self):
        r1 = noise_model({0:95, 1: 5})
        (n11,p11) = r1[1]
        r2 = noise_model({0:95, 1: 5, 2:0})
        (n21,p21) = r2[1]
        self.assertEqual(p11,p21)

    def test_large_numbers_at_the_end_make_earlier_numbers_more_probable(self):
        r1 = noise_model({0:95, 1: 5})
        (n11,p11) = r1[1]
        r2 = noise_model({0:95, 1: 5, 2:10})
        (n21,p21) = r2[1]
        self.assertGreater(p21, p11)

    def test_large_number_at_the_end_is_improbable(self):
        r = noise_model({0:75, 1: 5, 2:10})
        (n,p) = r[2]
        self.assertLess(p, 0.05)

    def test_small_number_at_the_end_is_probable(self):
        r = noise_model({0:944, 1: 50, 2:3})
        (n,p) = r[2]
        self.assertGreater(p, 0.05)

    def test_medium_number_at_the_end_is_vaguely_improbable(self):
        r = noise_model({0:941, 1: 50, 2:3, 3:1})
        (n,p) = r[3]
        self.assertLess(p, 0.05)
        self.assertGreater(p, 0.01)

if __name__ == '__main__':
    unittest.main()
