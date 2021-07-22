import unittest

import re
from marker_alignments.refdb_pattern import taxon_and_marker_patterns

def match_expected(self, refdb_format, input_string, expected_taxon, expected_marker):
    (taxon_p, marker_p) = taxon_and_marker_patterns(refdb_format)
    pattern_taxon = re.compile(taxon_p)
    taxon_search = pattern_taxon.search(input_string)
    taxon = next(g for g in taxon_search.groups() if g is not None)
    pattern_marker = re.compile(marker_p)
    marker_search = pattern_marker.search(input_string)
    marker = next(g for g in marker_search.groups() if g is not None)
    self.assertEqual(taxon, expected_taxon, str(taxon_search.groups()))
    self.assertEqual(marker, expected_marker, str( marker_search.groups()))


class ParseArgs(unittest.TestCase):

    def test_no_split(self):
        match_expected(self, "no-split", "", "", "")
        match_expected(self, "no-split", "xyz", "", "xyz")

    def test_eukprot(self):
        match_expected(self, "eukprot", "protist-Piridium_sociabile-418107at2759-S1", "Piridium_sociabile", "418107at2759-S1")

    def test_chocophlan(self):
        match_expected(self, "chocophlan", "39777__C4FSF9__HMPREF9321_0278|k__Bacteria.p__Firmicutes.c__Negativicutes.o__Veillonellales.f__Veillonellaceae.g__Veillonella.s__Veillonella_atypica|UniRef90_C4FSF9|UniRef50_D6KRB8|993", "Veillonella_atypica", "UniRef90_C4FSF9")

    def test_generic(self):
        match_expected(self, "generic", "", "", "")
        match_expected(self, "generic", "xyz", "", "xyz")
        match_expected(self, "generic", "protist-Piridium_sociabile-418107at2759-S1", "Piridium_sociabile", "418107at2759-S1")
        # we would prefer: Veillonella_atypica 
        match_expected(self, "generic", "39777__C4FSF9__HMPREF9321_0278|k__Bacteria.p__Firmicutes.c__Negativicutes.o__Veillonellales.f__Veillonellaceae.g__Veillonella.s__Veillonella_atypica|UniRef90_C4FSF9|UniRef50_D6KRB8|993", "39777__C4FSF9__HMPREF9321_0278", "UniRef90_C4FSF9")



if __name__ == '__main__':
    unittest.main()
