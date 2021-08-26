from scipy.stats import binom

def cutoff_fit_for_noise_model(taxon_counts_with_num_markers, noise_threshold):
    model = fit_noise_model(taxon_counts_with_num_markers)
    ns_above_threshold = [(n,p) for (n,p) in model if p <= noise_threshold]
    if not ns_above_threshold:
        return None
    else:
        return ns_above_threshold[0]

# When running an alignment and treating each match as presence of a marker,
# there are always some false positives.
# This is why presence of multiple markers is needed to identify a taxon.
# 
# Unfortunately with enough false positive markers, we start to hit false positive taxa
# and the threshold of how many markers we require to identfy a taxon starts to go up.
#
# We address this by fitting a null model to the data, and returning the first value for which
# the null model no longer fits.
def fit_noise_model(taxon_counts_with_num_markers):
    total_num_markers = sum([j * taxon_counts_with_num_markers[j] for j in taxon_counts_with_num_markers])
    num_taxa = sum([taxon_counts_with_num_markers[j] for j in taxon_counts_with_num_markers])
    if num_taxa == 0: return []

    # suppose (somewhat pessimistically) that there is no information content between markers found
    # and taxa present - the data would be just as good if each marker was independently assigned to a taxon
    # uniformly at random.
    # Then the number of markers for each taxon can be modelled as a binomially distributed random variable B_i, with:
    markers_for_each_taxon_n = total_num_markers
    markers_for_each_taxon_p = 1.0 / num_taxa

    results = []
    for num_markers in sorted(taxon_counts_with_num_markers.keys()):
        # there were this many taxa in the dataset
        taxon_count = taxon_counts_with_num_markers[num_markers]
        # null hypothesis: this is not yet above what we expect randomly
        # what's the chance of this happening?
        p = probability_at_least_taxon_count_num_markers_taxa(
                markers_for_each_taxon_n, markers_for_each_taxon_p,
                num_taxa, num_markers, taxon_count)
        results.append((num_markers, p))
    return results

def probability_at_least_taxon_count_num_markers_taxa(markers_for_each_taxon_n, markers_for_each_taxon_p, num_taxa, num_markers, taxon_count):
    if num_markers == 0:
        return 1
    if taxon_count == 0:
        return 1
    if num_taxa == 0 and num_markers > 0:
        return 0


    # calculate P(B_i = num_markers)
    num_markers_pmf = binom.pmf(k = num_markers, n = markers_for_each_taxon_n, p = markers_for_each_taxon_p)
    # define I(i, num_markers) to be 1 if B_i = num_markers and 0 otherwise
    # and let C(num_markers) be the sum of I(i, num_markers) over i
    # if I(i, num_markers) were independent, each C(num_markers) would follow a binomial distribution, with:
    taxa_with_num_markers_n = num_taxa
    taxa_with_num_markers_p = num_markers_pmf
    # they are not independent because sum of B_i adds up to total_num_markers but should be fine for small num_markers!

    # calculate the result using the survival function from the package:
    # P(C >= x) = 1 - cdf(x-1) = sf(x-1)
    sf = binom.sf(k = taxon_count - 1, n = taxa_with_num_markers_n, p = taxa_with_num_markers_p)
    return sf

