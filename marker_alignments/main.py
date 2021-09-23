import argparse
import sys
import pysam
import re
import math

from marker_alignments.store import AlignmentStore
from marker_alignments.refdb_pattern import taxon_and_marker_patterns

from marker_alignments.pysam2 import compute_contribution_to_marker_coverage, compute_alignment_identity


def next_g(search):
    return next(g for g in search.groups() if g is not None);

def taxon_and_marker(reference_name, pattern_taxon, pattern_marker, marker_to_taxon):

    taxon = marker_to_taxon[reference_name] if reference_name in marker_to_taxon else None
    taxon_search = pattern_taxon.search(reference_name)
    if taxon_search and taxon:
        taxon = "|".join([taxon, next_g(taxon_search)])
    elif taxon_search:
        taxon = next_g(taxon_search)
    elif taxon:
        taxon = taxon
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

def read_alignments(alignment_file, sqlite_db_path, pattern_taxon, pattern_marker, marker_to_taxon, min_mapq, min_query_length, min_match_identity):

    alignment_store = AlignmentStore(db_path=sqlite_db_path)
    alignment_store.start_bulk_write()

    for read in alignment_file.fetch():
        identity = compute_alignment_identity(read)
        if read.mapq < min_mapq:
            continue
        if read.infer_query_length() < min_query_length:
            continue
        if identity < min_match_identity:
            continue
        if not read.reference_name:
            raise ValueError("Read missing reference name: " + str(read))
        (taxon, marker) = taxon_and_marker(read.reference_name, pattern_taxon, pattern_marker, marker_to_taxon)
        if not taxon:
            raise ValueError("Could not find taxon in reference name: " + read.reference_name)
        if not marker:
            raise ValueError("Could not find marker in reference name: " + read.reference_name)

        alignment_store.add_alignment(
          taxon = taxon,
          marker = marker,
          query = read.query_name,
          identity = identity,
          coverage = compute_contribution_to_marker_coverage(alignment_file, read),
        )

    alignment_store.end_bulk_write()
    return alignment_store

def read_marker_to_taxon(path):
    result = {}
    with open(path, 'r') as f:
        for line in f:
            (marker, taxon) = line.rstrip().split("\t")
            result[marker] = taxon
    return result

output_type_options = ["marker_coverage", "marker_read_count", "marker_cpm", "marker_all", "taxon_coverage", "taxon_read_and_marker_count", "taxon_cpm", "taxon_all"]

def get_output(alignment_store, output_type, num_reads):
    if output_type == "marker_coverage":
        header = ["taxon", "marker", "marker_coverage"]
        lines = alignment_store.as_marker_coverage()
    elif output_type == "marker_read_count":
        header = ["taxon", "marker", "marker_read_count", "marker_avg_identity"]
        lines = alignment_store.as_marker_read_count()
    elif output_type == "marker_cpm":
        header = ["taxon", "marker", "marker_cpm"]
        lines = alignment_store.as_marker_cpm(num_reads)
    elif output_type == "marker_all":
        header = ["taxon", "marker", "marker_coverage", "marker_cpm", "marker_read_count", "marker_avg_identity"]
        lines = alignment_store.as_marker_all(num_reads)
    elif output_type == "taxon_coverage":
        header = ["taxon", "coverage"]
        lines = alignment_store.as_taxon_coverage()
    elif output_type == "taxon_read_and_marker_count":
        header = ["taxon", "taxon_num_reads", "taxon_num_markers", "taxon_max_reads_in_marker"]
        lines = alignment_store.as_taxon_read_and_marker_count()
    elif output_type == "taxon_cpm":
        header = ["taxon", "cpm"]
        lines = alignment_store.as_taxon_cpm(num_reads)
    elif output_type == "taxon_all":
        header = ["taxon", "coverage", "cpm", "taxon_num_reads", "taxon_num_markers", "taxon_max_reads_in_marker"]
        lines = alignment_store.as_taxon_all(num_reads)
    else:
        raise ValueError("Unknown output type: " + output_type)

    return header, lines

def main(argv=sys.argv[1:]):
    parser = argparse.ArgumentParser(
      description="summarize_marker_alignments - process and summarise alignments of metagenomic sequencing reads to reference databases of marker genes",
      formatter_class = argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--input", type=str, action="store", dest="input_alignment_file", help = "Input SAM/BAM", required=True)
    parser.add_argument("--sqlite-db-path", type=str, action="store", dest="sqlite_db_path", help = "Store a sqlite database under this path instead of in memory", default=None)
    parser.add_argument("--refdb-format", type=str, action="store", dest="refdb_format", help = "Reference database used for alignment, required for parsing reference names. Supported values: eukprot, chocophlan, generic, no-split (no split into marker and taxon)", default="generic")
    parser.add_argument("--refdb-regex-taxon", type=str, action="store", dest="refdb_regex_taxon", help = "Regex to read taxon name from reference name")
    parser.add_argument("--refdb-regex-marker", type=str, action="store", dest="refdb_regex_marker", help = "Regex to read marker name from reference name")
    parser.add_argument("--refdb-marker-to-taxon-path", type=str, action="store", dest="refdb_marker_to_taxon_path", help = "Lookup file, two columns - marker name, taxon name")
    parser.add_argument("--num-reads", type=int, action="store", dest="num_reads", help = "Total number of reads (required for CPM output)")
    parser.add_argument("--output-type", type=str, action="store", dest="output_type", help = "output type: "+", ".join(output_type_options), default = "marker_coverage")
    parser.add_argument("--output", type=str, action="store", dest="output_path", help = "output path", required=True)
    parser.add_argument("--min-read-mapq", type=int, action="store", dest="min_read_mapq", help = "when reading the input, skip alignments with MAPQ < min-read-mapq", default=0)
    parser.add_argument("--min-read-query-length", type=int, action="store", dest="min_read_query_length", help = "when reading the input, skip alignments shorter than min-read-query-length", default=0)
    parser.add_argument("--min-read-match-identity", type=float, action="store", dest="min_read_match_identity", help = "when reading the input, skip alignments where the proportion of matching bases in the alignment is less than min-read-match-identity", default=0)
    parser.add_argument("--min-taxon-num-markers", type=int, action="store", dest="min_taxon_num_markers", help = "For taxon output: only report taxa with at least min-taxon-num-markers markers")
    parser.add_argument("--min-taxon-num-reads", type=int, action="store", dest="min_taxon_num_reads", help = "For taxon output: only report taxa with at least min-taxon-num-reads reads")
    parser.add_argument("--min-taxon-fraction-primary-matches", type=float, action="store", dest="min_taxon_fraction_primary_matches", help = "Only keep taxa where no more than min-taxon-fraction-primary-matches fraction of alignments is inferior / secondary")
    parser.add_argument("--min-taxon-avg-match-identity", type=float, action="store", dest="min_taxon_avg_identity", help = "Only keep taxa where the average read identity is at least min-taxon-avg-match-identity")

    options=parser.parse_args(argv)

    if options.min_read_mapq and options.min_taxon_fraction_primary_matches:
        raise ValueError("It us unwise to combine --min-read-mapq and --min-taxon-fraction-primary-matches filters")

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

    if options.min_taxon_num_markers and options.output_type not in ["taxon_read_and_marker_count", "taxon_all"]:
        raise ValueError("--min-taxon-num-markers only valid for output types: taxon_read_and_marker_count, taxon_all")

    if options.min_taxon_num_reads and options.output_type not in ["taxon_read_and_marker_count", "taxon_all"]:
        raise ValueError("--min-taxon-num-reads only valid for output types: taxon_read_and_marker_count, taxon_all")

    alignment_store = read_alignments(
      alignment_file = pysam.AlignmentFile(options.input_alignment_file),
      sqlite_db_path = options.sqlite_db_path,
      marker_to_taxon = read_marker_to_taxon(options.refdb_marker_to_taxon_path) if options.refdb_marker_to_taxon_path else {},
      pattern_taxon = re.compile(options.refdb_regex_taxon),
      pattern_marker = re.compile(options.refdb_regex_marker),
      min_mapq = options.min_read_mapq,
      min_query_length = options.min_read_query_length,
      min_match_identity = options.min_read_match_identity,
    )

    if options.min_taxon_fraction_primary_matches:
        alignment_store.modify_table_filter_taxa_on_multiple_matches(min_fraction_primary_matches = options.min_taxon_fraction_primary_matches)

    if options.min_taxon_num_markers or options.min_taxon_num_reads:
        alignment_store.modify_table_filter_taxa_on_num_markers_and_reads(min_num_markers = options.min_taxon_num_markers or 0, min_num_reads = options.min_taxon_num_reads or 0)

    if options.min_taxon_avg_identity:
        alignment_store.modify_table_filter_taxa_on_avg_identity(min_avg_identity = options.min_taxon_avg_identity)

    header, lines = get_output(alignment_store, options.output_type, options.num_reads)

    field_formats = {
      "taxon" : "",
      "marker": "",
      "marker_cpm": ":.6f",
      "marker_coverage": ":.6f",
      "marker_read_count": ":.2f",
      "marker_avg_identity": ":.6f",
      "cpm" : ":.6f",
      "coverage" : ":.6f",
      "taxon_num_reads": ":.6f",
      "taxon_num_markers": ":d",
      "taxon_max_reads_in_marker": ":.6f",
    }
    formatter="\t".join(['{' + field_formats[field] +'}' for field in header]) + "\n"
    with open(options.output_path, 'w') as f:
        f.write("\t".join(header) + "\n")
        for line in lines:
            f.write(formatter.format(*line))

