import unittest

import logging
#logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('tests/filter_noise_model')
from marker_alignments.filter.noise_model import fit_noise_model, log_likelihood, cutoff_fit_for_noise_model, counts_as_list

class CountsAsList(unittest.TestCase):
    def test_cases(self):
        self.assertEqual(counts_as_list({0:1000}, 1), [1000, 0])
        self.assertEqual(counts_as_list({0:1000, 1: 100}, 2), [1000, 100, 0])
        self.assertEqual(counts_as_list({0:1000, 1: 100}, 1), [1100, 0, 0])
        self.assertEqual(counts_as_list({0:1000, 10: 100}, 1, length_limit = 3), [1100, 0, 0, 0])
        self.assertEqual(counts_as_list({0:1000, 10: 100}, 20, length_limit = 3), [1000, 0, 0, 100])


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

if __name__ == '__main__':
    unittest.main()
