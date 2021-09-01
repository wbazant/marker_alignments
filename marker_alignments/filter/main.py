import argparse
import sys
import csv

import logging
logging.basicConfig(format='%(asctime)s filter_marker_alignments_taxa: %(message)s', datefmt='%d-%b-%y %H:%M:%S')
logger=logging.getLogger(__name__)

from marker_alignments.filter.noise_model import cutoff_fit_for_noise_model

def main(argv=sys.argv[1:]):
    parser = argparse.ArgumentParser(
      description="filter_marker_alignments_taxa - apply a filter to taxon output of summarize_marker_alignments",
      formatter_class = argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--input", type=str, action="store", dest="input_path", help = "Input summary file", required=True)
    parser.add_argument("--require-min-markers", type=int, action="store", dest="require_min_markers", help = "Require min markers to keep a taxon")
    parser.add_argument("--use-noise-model-for-min-markers", action="store_true", dest="use_noise_model_for_min_markers", help = "Use a null model where markers associate with taxa at random, and select the most appropriate value for --require-min-markers")
    parser.add_argument("--total-num-taxa", type=int, action="store", dest="total_num_taxa", help = "Total number of taxa in the reference - required for fitting the noise model")
    parser.add_argument("--taxon-to-markers-beta-sample-size", type=float, action="store", dest="beta_sample_size", help = "Sample size (sum of shape parameters a and b when proportion of markers per taxon is modelled as a beta distribution) - required for fitting the noise model")
    parser.add_argument("--output", type=str, action="store", dest="output_path", help = "output path", required=True)
    parser.add_argument("--verbose", action="store_true", dest="verbose", help = "turn on logging")

    options=parser.parse_args(argv)
    if options.verbose:
        logger.setLevel(logging.INFO)

    if options.use_noise_model_for_min_markers and not options.total_num_taxa:
        raise ValueError("--total-num-taxa required if --use-noise-model-with-confidence-threshold is provided")

    if options.use_noise_model_for_min_markers and not options.beta_sample_size:
        raise ValueError("--taxon-to-markers-beta-sample-size required if --use-noise-model-with-confidence-threshold is provided")

    lines, fieldnames = [], {}
    with open(options.input_path) as f:
        reader = csv.DictReader(f, delimiter='\t', quoting=csv.QUOTE_NONE)
        lines = [ l for l in reader ]
        fieldnames = reader.fieldnames
    
    logger.info("tab-separated file with %s data rows read from %s", len(lines), options.input_path)

    if options.total_num_taxa and options.total_num_taxa < len(lines):
        raise ValueError("--total-num-taxa provided (%s) is lower than the number of data rows (%s)".format(options.total_num_taxa, len(lines)))

    require_min_markers = None
    if options.use_noise_model_for_min_markers:
        taxon_counts_with_num_markers = {0 : options.total_num_taxa - len(lines)}
        for l in lines:
            n = int(l["taxon_num_markers"])
            if n not in taxon_counts_with_num_markers:
                taxon_counts_with_num_markers[n] = 0
            taxon_counts_with_num_markers[n] += 1
        cutoff_fit = cutoff_fit_for_noise_model(taxon_counts_with_num_markers, options.beta_sample_size, logger)
        if options.require_min_markers and options.require_min_markers > cutoff_fit:
            logger.info("Cutoff fit is less than --require-min-markers value - will use that instead: %s", options.require_min_markers)
            require_min_markers = options.require_min_markers
        else:
            require_min_markers = cutoff_fit
    elif options.require_min_markers:
        require_min_markers = options.require_min_markers

    if require_min_markers:
        lines = [ l for l in lines if int(l["taxon_num_markers"]) >= require_min_markers ]
        logger.info("Kept %s taxa with at least %s markers", len(lines), require_min_markers)


    with open(options.output_path, 'w') as f:
        writer = csv.DictWriter(f, delimiter='\t', fieldnames = fieldnames)
        writer.writeheader()
        for l in lines:
            writer.writerow(l)

    logger.info("tab-separated file with %s data rows written to %s", len(lines), options.output_path)


