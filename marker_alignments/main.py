import argparse
import sys
import pysam
import re
import math

from marker_alignments.store import AlignmentStore
from marker_alignments.write import write, output_type_options
from marker_alignments.mcl import clusters 
from marker_alignments.refdb_pattern import taxon_and_marker_patterns

from marker_alignments.pysam2 import compute_contribution_to_marker_coverage, compute_alignment_identity

def next_g(search):
    return next(g for g in search.groups() if g is not None)

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

def parse_arguments(argv=sys.argv[1:]):
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
    parser.add_argument("--min-taxon-num-markers", type=int, action="store", dest="min_taxon_num_markers", help = "Only keep taxa with at least min-taxon-num-markers markers")
    parser.add_argument("--min-taxon-num-reads", type=int, action="store", dest="min_taxon_num_reads", help = "Only keep taxa with at least min-taxon-num-reads reads")
    parser.add_argument("--min-taxon-fraction-primary-matches", type=float, action="store", dest="min_taxon_fraction_primary_matches", help = "Only keep taxa where no more than min-taxon-fraction-primary-matches fraction of alignments is inferior / secondary")
    parser.add_argument("--min-taxon-better-marker-cluster-averages-ratio", type=float, action="store", dest="min_taxon_better_cluster_averages_ratio", help = "Only keep taxa where the ratio between markers which have at least average match identity relative to their clusters and markers with identity below average is at least min-taxon-better-cluster-averages-ratio")

    parser.add_argument("--threshold-avg-match-identity-to-call-known-taxon", type=float, action="store", dest="threshold_identity_to_call_taxon", help = "Threshold on average match identity to return taxon in reference")
    parser.add_argument("--threshold-num-reads-to-call-unknown-taxon", type=int, action="store", dest="threshold_num_reads_to_call_unknown_taxon", help = "To positively identify an unknown taxon (fits all criteria except match identity) expect this many reads from a taxon cluster")
    parser.add_argument("--threshold-num-markers-to-call-unknown-taxon", type=int, action="store", dest="threshold_num_markers_to_call_unknown_taxon", help = "To positively identify an unknown taxon (fits all criteria except match identity) expect this many markers from a taxon cluster")
    parser.add_argument("--threshold-num-taxa-to-call-unknown-taxon", type=int, action="store", dest="threshold_num_taxa_to_call_unknown_taxon", help = "To positively identify an unknown taxon (fits all criteria except match identity) expect this many taxa from a taxon cluster")
    result = parser.parse_args(argv)
    return result

def main(argv=sys.argv[1:]):
    options=parse_arguments(argv)

    if options.min_read_mapq and (options.min_taxon_fraction_primary_matches or options.min_taxon_better_cluster_averages_ratio) :
        raise ValueError("It us unwise to combine --min-read-mapq and filters that rely on secondary matches!")

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
      marker_to_taxon = read_marker_to_taxon(options.refdb_marker_to_taxon_path) if options.refdb_marker_to_taxon_path else {},
      pattern_taxon = re.compile(options.refdb_regex_taxon),
      pattern_marker = re.compile(options.refdb_regex_marker),
      min_mapq = options.min_read_mapq,
      min_query_length = options.min_read_query_length,
      min_match_identity = options.min_read_match_identity,
    )

    alignment_store.cluster_markers_by_matches()

    if options.min_taxon_better_cluster_averages_ratio:
        alignment_store.modify_table_filter_taxa_on_cluster_averages(min_better_cluster_averages_ratio = options.min_taxon_better_cluster_averages_ratio)

    if options.min_taxon_fraction_primary_matches:
        alignment_store.modify_table_filter_taxa_on_multiple_matches(min_fraction_primary_matches = options.min_taxon_fraction_primary_matches)

    if options.min_taxon_num_markers or options.min_taxon_num_reads:
        alignment_store.modify_table_filter_taxa_on_num_markers_and_reads(min_num_markers = options.min_taxon_num_markers or 0, min_num_reads = options.min_taxon_num_reads or 0)


    alignment_store.cluster_taxa_by_matches()

    if options.threshold_identity_to_call_taxon or options.threshold_num_reads_to_call_unknown_taxon or options.threshold_num_markers_to_call_unknown_taxon or options.threshold_num_taxa_to_call_unknown_taxon:
        alignment_store.modify_table_transform_taxa_on_thresholds_and_clusters(
                threshold_identity = options.threshold_identity_to_call_taxon or 0,
                min_num_taxa_below_identity = options.threshold_num_taxa_to_call_unknown_taxon or 0,
                min_num_markers_below_identity = options.threshold_num_markers_to_call_unknown_taxon or 0,
                min_num_reads_below_identity = options.threshold_num_reads_to_call_unknown_taxon or 0
        )

    write(alignment_store, options.output_type, options.output_path, options.num_reads)
