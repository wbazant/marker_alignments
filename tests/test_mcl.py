import unittest
import pysam

from marker_alignments.mcl import clusters


class Mcl(unittest.TestCase):

    def test_null_case(self):
        self.assertEqual(clusters([]), [])

    def test_close_tie(self):
        result = clusters([
          ["t1", "t2", 10],
          ["t3", "t4", 10],
          ["t3", "t5", 10],
          ["t1", "t5", 10],
          ["t2", "t5", 10],
          ["t3", "t5", 10],
          ["t4", "t5", 10.1],
        ])
        self.assertEqual(result, [['t3', 't4', 't5'], ['t1', 't2']])

        result = clusters([
          ["t1", "t2", 10],
          ["t3", "t4", 10],
          ["t3", "t5", 10],
          ["t1", "t5", 10],
          ["t2", "t5", 10],
          ["t3", "t5", 10],
          ["t4", "t5", 9.9],
        ])
        self.assertEqual(result,[['t1', 't2', 't5'], ['t3', 't4']])

if __name__ == '__main__':
    unittest.main()
