import argparse
import sys
import statistics

from scipy.stats import beta

def read_taxon_to_num_markers(path):
    taxon_to_num_markers = {}
    with open(path, 'r') as f:
        for line in f:
            (marker, taxon) = line.rstrip().split("\t")
            if taxon not in taxon_to_num_markers:
                taxon_to_num_markers[taxon] = 0
            taxon_to_num_markers[taxon]+=1

    return taxon_to_num_markers

def main(argv=sys.argv[1:]):
    parser = argparse.ArgumentParser(
      description="refdb stats",
      formatter_class = argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--input", type=str, action="store", dest="input_refdb_path", help = "Input refdb", required=True)

    options=parser.parse_args(argv)

    taxon_to_num_markers = read_taxon_to_num_markers(options.input_refdb_path)

    print("Num taxa: {}".format(len(taxon_to_num_markers.keys())))
    total_num_markers = sum(taxon_to_num_markers.values())
    print("Num markers: {}".format(total_num_markers))
    print("Mean markers: {}".format(statistics.mean(taxon_to_num_markers.values())))
    print("Variance markers: {}".format(statistics.variance(taxon_to_num_markers.values())))

    num_markers_fracs = [1.0 * x / total_num_markers for x in sorted(taxon_to_num_markers.values())]
    (a, b, loc, scale) = beta.fit(num_markers_fracs, floc=0, fscale=1)
    beta_mean = a /(a+b)
    beta_mode = (a-1)/(a+b -2)
    beta_sample_size = a+b
    print("Shape parameters for a beta fit:\n  a = {}\n  b = {}\n  mean = {}\n  mode = {}\n  sample size = {}".format(a,b, beta_mean, beta_mode, beta_sample_size))

if __name__ == "__main__":
    main()
