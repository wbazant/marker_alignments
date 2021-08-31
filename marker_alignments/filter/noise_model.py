import numpy
from scipy.stats import binom, betabinom, multinomial

def cutoff_fit_for_noise_model(taxon_counts_with_num_markers, beta_sample_size, logger):
    m = max(taxon_counts_with_num_markers.keys())
    if m < 2:
        return 2

    log_likelihoods = []
    for candidate_cutoff in candidate_cutoffs(taxon_counts_with_num_markers):
        ks = counts_as_list(taxon_counts_with_num_markers, candidate_cutoff)

        ll = log_likelihood(ks, beta_sample_size)
        logger.info("Cutoff %s: log likelihood %s", candidate_cutoff, ll)
        log_likelihoods.append((ll, -candidate_cutoff))

        fit_noise_model({j: ks[j] for j in range(0, min(20,m+2) if m < 20 else 21)}, beta_sample_size, logger)

    logger.info("Likelihoods of different possible cutoffs:\nCutoff\tLog likelihood\n%s", "\n".join(["{}\t{:6f}".format(-cutoff, ll) for (ll, cutoff) in log_likelihoods]))
    log_likelihoods.sort(reverse=True)
    (lll, cutoff) = log_likelihoods[0]
    return -cutoff

def counts_as_list(taxon_counts_with_num_markers, candidate_cutoff, length_limit = 20):
    m = max(taxon_counts_with_num_markers.keys())
    ks = []
    for j in range(0, min(length_limit, m+2)):
       if j in taxon_counts_with_num_markers and j < candidate_cutoff:
           ks.append(taxon_counts_with_num_markers[j])
       else:
           if j in taxon_counts_with_num_markers:
               ks[0]+=taxon_counts_with_num_markers[j]
           ks.append(0)

    if m >= length_limit:
        k_last = 0
        for jj in range(length_limit, m + 2):
           if jj in taxon_counts_with_num_markers and jj < candidate_cutoff:
               k_last += taxon_counts_with_num_markers[jj]
           elif jj in taxon_counts_with_num_markers:
               ks[0]+=taxon_counts_with_num_markers[jj]
        ks.append(k_last)
    return ks

# for each value in the dataset, consider making that the last value below cutoff
def candidate_cutoffs(taxon_counts_with_num_markers):
    result = set()
    for x in taxon_counts_with_num_markers:
        if x > 0:
            result.add(x+1)
    return sorted(result)

# Suppose the number of markers for each taxon followed a beta binomial distribution
# with sample size as provided, and mean set to the average number of markers
# num_taxa independent trials, questions about the sum - that's a multinomial distribution
def log_likelihood(ks, beta_sample_size):
    total_num_markers = sum([j * ks[j] for j in range(0, len(ks))])
    num_taxa = sum(ks)

    markers_for_each_taxon_n = total_num_markers
    markers_for_each_taxon_p = 1.0 / num_taxa

    markers_for_each_taxon_shape_a = markers_for_each_taxon_p * beta_sample_size
    markers_for_each_taxon_shape_b = (1 - markers_for_each_taxon_p) * beta_sample_size

    markers_for_each_taxon_rv = betabinom(n = markers_for_each_taxon_n, a = markers_for_each_taxon_shape_a, b = markers_for_each_taxon_shape_b)

    ps = [markers_for_each_taxon_rv.pmf(k = k) for k in range(0, len(ks))] 

    ks.reverse()
    ps.reverse()
    ll = multinomial.logpmf(x = ks, n=num_taxa, p = ps)
    if numpy.isnan(ll):
        ll = numpy.NINF
    ks.reverse()
    ps.reverse()
    return ll


# When running an alignment and treating each match as presence of a marker,
# there are always some false positives.
# This is why presence of multiple markers is needed to identify a taxon.
# 
# Unfortunately with enough false positive markers, we start to hit false positive taxa
# and the threshold of how many markers we require to identfy a taxon starts to go up.
def fit_noise_model(taxon_counts_with_num_markers, beta_sample_size, logger):
    total_num_markers = sum([j * taxon_counts_with_num_markers[j] for j in taxon_counts_with_num_markers])
    num_taxa = sum([taxon_counts_with_num_markers[j] for j in taxon_counts_with_num_markers])
    if num_taxa == 0: return []

    # suppose (somewhat pessimistically) that there is no information content between markers found
    # and taxa present - the data would be just as good if each marker was independently assigned to a taxon
    # uniformly at random.
    # Then the number of markers for each taxon B_i can be modelled as a binomially distributed random variable, with:
    markers_for_each_taxon_n = total_num_markers
    markers_for_each_taxon_p = 1.0 / num_taxa

    # That actually doesn't fit super well because taxa in the reference have different sizes - 
    # if we pick a subset of markers at random, some taxa will get more of them than others.
    # Instead, model B_i as a beta binomial distributed random variable with the same mean.

    # https://en.wikipedia.org/wiki/Beta_distribution#Mean_and_sample_size
    markers_for_each_taxon_shape_a = markers_for_each_taxon_p * beta_sample_size
    markers_for_each_taxon_shape_b = (1 - markers_for_each_taxon_p) * beta_sample_size

    log=["Distributing {} markers across {} taxa, counts of taxa with k markers are:".format(total_num_markers, num_taxa)]
    log.append("\t".join(["k", "actual", "expected", "pmf", "p(at least actual)"]))
    results = []
    for num_markers in sorted(taxon_counts_with_num_markers.keys()):
        # there were this many taxa in the dataset
        taxon_count = taxon_counts_with_num_markers[num_markers]
        # null hypothesis: this is not yet above what we expect randomly
        # what's the chance of this happening?

        num_markers_pmf = betabinom.pmf(k = num_markers, n = markers_for_each_taxon_n, a = markers_for_each_taxon_shape_a, b = markers_for_each_taxon_shape_b)
        # num_markers_pmf = binom.pmf(k = num_markers, n = markers_for_each_taxon_n, p = markers_for_each_taxon_p)
        p = probability_at_least_taxon_count_num_markers_taxa(
                num_markers_pmf,
                num_taxa, num_markers, taxon_count)
        log.append("%s%s\t%s\t%.2f\t%.2g\t%.2g" % (num_markers, "+" if num_markers == max(taxon_counts_with_num_markers.keys()) else "", taxon_count, round(num_markers_pmf * num_taxa, 2), num_markers_pmf, p ))
        results.append((num_markers, p))
    logger.info("\n".join(log))
    return results

def probability_at_least_taxon_count_num_markers_taxa(num_markers_pmf, num_taxa, num_markers, taxon_count):
    if num_markers == 0:
        return 1
    if taxon_count == 0:
        return 1
    if num_taxa == 0 and num_markers > 0:
        return 0

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

