
eukprot_refdb_regex_taxon="^[^-]+-([^-]+)-.*$"
eukprot_refdb_regex_marker="^[^-]+-[^-]+-(.*)$"

chocophlan_refdb_regex_taxon="s__(.*?)\\|"
chocophlan_refdb_regex_marker="(UniRef90[^|]*)"


no_split_regex_taxon="(^)"
no_split_regex_marker="(.*)"

def taxon_and_marker_patterns(refdb_format):
    if refdb_format == "eukprot":
        return (eukprot_refdb_regex_taxon, eukprot_refdb_regex_marker)
    elif refdb_format == "chocophlan":
        return (chocophlan_refdb_regex_taxon, chocophlan_refdb_regex_marker)
    elif refdb_format == "no-split":
        return (no_split_regex_taxon, no_split_regex_marker)
    elif refdb_format == "generic":
        return (
           "|".join([eukprot_refdb_regex_taxon, chocophlan_refdb_regex_taxon, "(^[^:|]*)[:|]", no_split_regex_taxon]),
           "|".join([eukprot_refdb_regex_marker, chocophlan_refdb_regex_marker, "[:|]?([^:|]*)$"]),
        )
    else:
        return (None, None)

