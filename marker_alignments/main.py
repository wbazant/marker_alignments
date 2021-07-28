import argparse
import sys
import pysam
import re
import math

from marker_alignments.store import AlignmentStore
from marker_alignments.refdb_pattern import taxon_and_marker_patterns

# if a read matches multiple markers, distribute it across matches
# do it proportionally to the second power of match identity
# (identity = match length / alignment length)
def multiple_matches_weighing(read): #pysam.AlignedSegment
    match_length = len(read.get_aligned_pairs(matches_only=True))
    alignment_length = len(read.get_aligned_pairs())
    identity_fraction = 1.0 * match_length / alignment_length
    return round(math.pow(identity_fraction, 2), 6)

# longer genes have more reads aligning to them
# so instead of counting reads, estimate a number of copies of each marker
# this should be proportional to species abundance and agree between markers
# ( and also proportional to sequencing depth - we correct that later if calculating CPM output)
def contribution_to_marker_coverage(read, marker_length):
    return round(1.0 * read.infer_query_length() / marker_length, 6)

def next_g(search):
    return next(g for g in search.groups() if g is not None);

def taxon_and_marker(reference_name, pattern_taxon, pattern_marker, marker_to_taxon_id):
    taxon_id = marker_to_taxon_id[reference_name] if reference_name in marker_to_taxon_id else None

    taxon_search = pattern_taxon.search(reference_name)
    if taxon_search and taxon_id:
        taxon = "|".join([taxon_id, next_g(taxon_search)])
    elif taxon_search:
        taxon = next_g(taxon_search)
    elif taxon_id:
        taxon = taxon_id
    else:
        taxon = None

    marker_search = pattern_marker.search(reference_name)
    if marker_search:
        marker=next_g(marker_search)
    elif taxon_search:
        # Regexes seem to work, just not this one 
        marker = None
    else:
        marker = reference_name
    return (taxon, marker)


def read_alignments(alignment_file, sqlite_db_path, pattern_taxon, pattern_marker, marker_to_taxon_id):

    alignment_store = AlignmentStore(db_path=sqlite_db_path)
    alignment_store.start_bulk_write()

    for read in alignment_file.fetch():
        (taxon, marker) = taxon_and_marker(read.reference_name, pattern_taxon, pattern_marker, marker_to_taxon_id)
        if not taxon:
            raise ValueError("Could not find taxon in reference name: " + read.reference_name)
        if not marker:
            raise ValueError("Could not find marker in reference name: " + read.reference_name)

        alignment_store.add_alignment(
          taxon = taxon,
          marker = marker,
          query = read.query_name,
          weight = multiple_matches_weighing(read),
          coverage = contribution_to_marker_coverage(read, alignment_file.get_reference_length(read.reference_name)),
        )

    alignment_store.end_bulk_write()
    return alignment_store

def read_marker_to_taxon_id(path):
    result = {}
    with open(path, 'r') as f:
        for line in f:
            (marker, taxon_id) = line.rstrip().split("\t")
            result[marker] = taxon_id
    return result

output_type_options = ["marker_coverage", "marker_read_count", "marker_cpm", "marker_all", "taxon_coverage", "taxon_read_and_marker_count", "taxon_cpm", "taxon_all"]

def main(argv=sys.argv[1:]):
    parser = argparse.ArgumentParser(
      description="marker_alignment - process and summarise alignments of metagenomic sequencing reads to reference databases of marker genes",
      formatter_class = argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--input", type=str, action="store", dest="input_alignment_file", help = "Input SAM/BAM", required=True)
    parser.add_argument("--sqlite-db-path", type=str, action="store", dest="sqlite_db_path", help = "Store a sqlite database under this path instead of in memory", default=None)
    parser.add_argument("--refdb-format", type=str, action="store", dest="refdb_format", help = "Reference database used for alignment, required for parsing reference names. Supported values: eukprot, chocophlan, generic, no-split (no split into marker and taxon)", default="generic")
    parser.add_argument("--refdb-regex-taxon", type=str, action="store", dest="refdb_regex_taxon", help = "Regex to read taxon name from reference name")
    parser.add_argument("--refdb-regex-marker", type=str, action="store", dest="refdb_regex_marker", help = "Regex to read marker name from reference name")
    parser.add_argument("--refdb-marker-to-taxon-id-path", type=str, action="store", dest="refdb_marker_to_taxon_id_path", help = "Lookup file, two columns - marker name, taxon name")
    parser.add_argument("--num-reads", type=int, action="store", dest="num_reads", help = "Total number of reads (required for CPM output)")
    parser.add_argument("--output-type", type=str, action="store", dest="output_type", help = "output type: "+", ".join(output_type_options), default = "marker_coverage")
    parser.add_argument("--output", type=str, action="store", dest="output_path", help = "output path", required=True)

    options=parser.parse_args(argv)

    if options.output_type not in output_type_options:
        raise ValueError("Unknown output type: " + options.output_type)

    if options.refdb_format:
        (tp, mp) = taxon_and_marker_patterns(options.refdb_format)
        if not (tp and mp):
            raise ValueError("Unknown refdb format: " + options.refdb_format)
        options.refdb_regex_taxon=tp
        options.refdb_regex_marker=mp

    if not (options.refdb_regex_taxon and options.refdb_regex_marker):
        raise ValueError("Please provide either refdb format, or taxon + marker regexes")


    if options.output_type in ["marker_all","marker_cpm", "taxon_all", "taxon_cpm"] and not options.num_reads:
        raise ValueError("--num-reads required for calculating " + options.output_type)

    alignment_store = read_alignments(
      alignment_file = pysam.AlignmentFile(options.input_alignment_file),
      sqlite_db_path = options.sqlite_db_path,
      marker_to_taxon_id = read_marker_to_taxon_id(options.refdb_marker_to_taxon_id_path) if options.refdb_marker_to_taxon_id_path else {},
      pattern_taxon = re.compile(options.refdb_regex_taxon),
      pattern_marker = re.compile(options.refdb_regex_marker),
    )

    if options.output_type == "marker_coverage":
        header = ["taxon", "marker", "marker_coverage"]
        lines = alignment_store.as_marker_coverage()
    elif options.output_type == "marker_read_count":
        header = ["taxon", "marker", "marker_read_count"]
        lines = alignment_store.as_marker_read_count()
    elif options.output_type == "marker_cpm":
        header = ["taxon", "marker", "marker_cpm"]
        lines = alignment_store.as_marker_cpm(options.num_reads)
    elif options.output_type == "marker_all":
        header = ["taxon", "marker", "marker_coverage", "marker_cpm", "marker_read_count"]
        lines = alignment_store.as_marker_all(options.num_reads)
    elif options.output_type == "taxon_coverage":
        header = ["taxon", "coverage"]
        lines = alignment_store.as_taxon_coverage()
    elif options.output_type == "taxon_read_and_marker_count":
        header = ["taxon", "taxon_num_reads", "taxon_num_markers"]
        lines = alignment_store.as_taxon_read_and_marker_count()
    elif options.output_type == "taxon_cpm":
        header = ["taxon", "cpm"]
        lines = alignment_store.as_taxon_cpm(options.num_reads)
    elif options.output_type == "taxon_all":
        header = ["taxon", "coverage", "cpm", "taxon_num_reads", "taxon_num_markers"]
        lines = alignment_store.as_taxon_all(options.num_reads)

    field_formats = {
      "taxon" : "",
      "marker": "",
      "marker_cpm": ":.6f",
      "marker_coverage": ":.6f",
      "marker_read_count": ":.2f",
      "cpm" : ":.6f",
      "coverage" : ":.6f",
      "taxon_num_reads": ":.6f",
      "taxon_num_markers": ":d",
    }
    formatter="\t".join(['{' + field_formats[field] +'}' for field in header]) + "\n"
    with open(options.output_path, 'w') as f:
        f.write("\t".join(header) + "\n")
        for line in lines:
            f.write(formatter.format(*line))

