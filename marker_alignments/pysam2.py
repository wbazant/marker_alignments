import pysam
import itertools

# longer genes have more reads aligning to them
# so instead of counting reads, estimate a number of copies of each marker
# this should be proportional to species abundance and agree between markers
# ( and also proportional to sequencing depth - we correct that later if calculating CPM output)
def compute_contribution_to_marker_coverage(alignment_file, sam_record):
    marker_length = alignment_file.get_reference_length(sam_record.reference_name)

    return round(1.0 * sam_record.infer_query_length() / marker_length, 6)

# adapted from https://github.com/mortazavilab/TALON/blob/master/src/talon/transcript_utils.py
def compute_alignment_identity(sam_record):
    """ This function computes what fraction of the read matches the reference
        genome."""

    read_ID = sam_record.query_name
    try:
        MD_tag = sam_record.get_tag('MD')
    except KeyError:
        raise ValueError("SAM transcript %s lacks an MD tag" % read_ID)

    total_bases = sam_record.infer_query_length()
    matches = 0.0
    ops, counts = splitMD(MD_tag)
    for op,ct in zip(ops, counts):
        if op == "M":
            matches += ct
        if op == "D":
            total_bases += ct

    return round(1.0 * matches/total_bases, 6)

def splitMD(MD):
        """ Takes MD tag and splits into two lists:
            one with capital letters (match operators), and one with
            the number of bases that each operation applies to. """

        operations = []

        # Split MD string where type changes.
        # Digits are separated from base changes.
        # Deletions (with ^) are captured together.
        counts = ["".join(x) for _, x in itertools.groupby(MD, key=str.isdigit)]

        # Get operations
        for i in range(0,len(counts)):
            curr = counts[i]
            try:
                counts[i] = int(curr)
                operations.append("M")
            except ValueError:
                # Handle deletion
                if curr.startswith("^"):
                    operations.append("D")
                    counts[i] = len(counts[i]) - 1
                else:
                    operations.append("X")
                    counts[i] = len(counts[i])

        return operations, counts
