import unittest
import pysam

import sys
import os 
dir_path = os.path.dirname(os.path.realpath(__file__))
input_path = dir_path + "/data/example.sam"

import re
import tempfile

from marker_alignments.summarize.main import main, output_type_options

def read_lines_and_remove(path):
    if not os.path.isfile(path):
        return []

    with open(path) as f:
        lines = f.readlines()
    os.remove(path)
    return lines

def run(self, optional_args=[]):
    output_path = tempfile.mktemp()
    main(["--input", input_path, "--output", output_path, *optional_args])
    lines = read_lines_and_remove(output_path)
    self.assertTrue(len(lines) > 1, msg=lines)

def fail(self, exception, args):
    with self.assertRaises(exception):
        sys.stderr = open(os.devnull, 'w')
        sys.stdout = open(os.devnull, 'w')
        main(args)
        sys.stderr = sys.__stderr__
        sys.stdout = sys.__stdout__

class ParseArgs(unittest.TestCase):

    def test_fail_bad_args(self):
        output_path = tempfile.mktemp()
        fail(self, SystemExit, [])
        fail(self, SystemExit, ["-h"])
        fail(self, SystemExit, ["--input", input_path])
        fail(self, SystemExit, ["--output", output_path])
        fail(self, ValueError, ["--input", input_path, "--output", output_path, "--refdb-format", "x"])
        fail(self, ValueError, ["--input", input_path, "--output", output_path, "--output-type", "x"])
        fail(self, ValueError, ["--input", input_path, "--output", output_path, "--output-type", "marker_all"])

    def test_no_optional_args(self):
        run(self, [])

    def test_output_types(self):
        for output_type in output_type_options:
            with self.subTest(output_type):
                run(self, ["--output-type", output_type, "--num-reads", "42"])

if __name__ == '__main__':
    unittest.main()
