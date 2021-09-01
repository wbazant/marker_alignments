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
        log_likelihoods.append((ll, -candidate_cutoff))

    logger.info("Likelihoods of different possible cutoffs:\nCutoff\tLog likelihood\n%s", "\n".join(["{}\t{:6f}".format(-cutoff, ll) for (ll, cutoff) in log_likelihoods]))

    log_likelihoods.sort(reverse=True)
    (lll, cutoff) = log_likelihoods[0]
    cutoff = -cutoff
    logger.info("Selected a cutoff: %s", cutoff)

    ks = counts_as_list(taxon_counts_with_num_markers, cutoff)
    actuals = counts_as_list(taxon_counts_with_num_markers, 10000000000)

    d = {j: ks[j] for j in range(0, len(ks))}
    expected = fit_noise_model(d, beta_sample_size)
    log=["Expected under this cutoff:"]
    log.append("\t".join(["n", "actual", "expected"]))
    for j in range(0, len(ks)):
        if j == cutoff:
            log.append("---\t---\t---\t")
        log.append("{}{}\t{}\t{:.2f}".format(j, "+" if j == len(ks)-1 else "", actuals[j],  expected[j]))
    logger.info("%s", "\n".join(log))
    return cutoff

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
    if taxon_counts_with_num_markers:
        # also try 2
        # if there are no 1s but there are higher values, this will be a cutoff with likelihood 0
        result.add(2)

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


def fit_noise_model(taxon_counts_with_num_markers, beta_sample_size):
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

    results = {}
    for num_markers in sorted(taxon_counts_with_num_markers.keys()):
        taxon_count = taxon_counts_with_num_markers[num_markers]
        num_markers_pmf = betabinom.pmf(k = num_markers, n = markers_for_each_taxon_n, a = markers_for_each_taxon_shape_a, b = markers_for_each_taxon_shape_b)
        expected = 1.0 * num_taxa * num_markers_pmf
        results[num_markers] = expected
    return results
