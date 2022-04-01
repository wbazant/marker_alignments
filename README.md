# marker alignments
This a Python package to process and summarise alignments of metagenomic sequencing reads to a reference database of marker genes. You can use it in conjunction with an aligner like `bowtie2` to produce an estimate of taxa present in a metagenomic sample.

The package was developed in the context of looking for eukaryotes - most of the facilities are for producing good guesses from small amounts of potentially unreliable information. There are read level filters, clustering facilities for making sense of multiple alignments per query, and a number of thresholds.

## Installation

To install via pip:
```
pip install marker_alignments
```

## Usage

### Introduction

Download a small example alignment file, and run `marker_alignments` with most basic options:

```
wget "https://raw.githubusercontent.com/wbazant/marker_alignments/main/tests/data/example.sam"

marker_alignments --input example.sam --output /dev/stdout
```

If the package installed correctly, you should see a coverage report for each reference in the alignments file. `marker_alignments --help` should show you all filtering options.


### Detecting eukaryotes

First download the EukDetect reference database following [EukDetect installation instructions](https://github.com/allind/EukDetect).

Then follow this example to download an example metagenomic file, run alignments to a reference database bundled with EukDetect, and obtain a profile using suitable filtering options:

```
REFDB_LOCATION="eukdb"

wget "ftp.sra.ebi.ac.uk/vol1/fastq/ERR274/009/ERR2749179/ERR2749179_1.fastq.gz"
wget "ftp.sra.ebi.ac.uk/vol1/fastq/ERR274/009/ERR2749179/ERR2749179_2.fastq.gz"
gunzip *gz
FASTQ_1="ERR2749179_1.fastq"
FASTQ_2="ERR2749179_2.fastq"

bowtie2 --omit-sec-seq --no-discordant --no-unal \
  -x $REFDB_LOCATION/ncbi_eukprot_met_arch_markers.fna \
  -k10,10
  -1 ERR2749179_1.fastq.gz \
  -2 ERR2749179_2.fastq.gz \
  -S ERR2749179.sam 

FILTERING_OPTS="--min-read-query-length 60 --min-taxon-num-markers 2 --min-taxon-num-reads 2 --min-taxon-better-marker-cluster-averages-ratio 1.01 --threshold-avg-match-identity-to-call-known-taxon 0.97  --threshold-num-taxa-to-call-unknown-taxon 1 --threshold-num-markers-to-call-unknown-taxon 4     --threshold-num-reads-to-call-unknown-taxon 8"

marker_alignments --input ERR2749179.sam --output ERR2749179.taxa.tsv \
  --refdb-format eukprot \
  --refdb-marker-to-taxon-path $REFDB_LOCATION/busco_taxid_link.txt \
  --output-type taxon_all \
  --num-reads $(grep -c '^@' $FASTQ_1) \
  $FILTERING_OPTS
```

To do this for multiple samples, try the Nextflow pipeline [wbazant/CORRAL](https://github.com/wbazant/CORRAL).


### Other uses
The basic workflow type supported by this package is to give it an alignment file, and look at reports produced. 

There are multiple filtering options aiming to reduce noise enough that the resulting taxonomic profile can be passed on to other tools. Alternatively, you can specify `--output-type pairs_of_taxa_shared_queries` or `output-type taxa_in_marker_clusters` and look at shared alignments between the queries, and get a detailed view of what the sequences in your metagenomic sample are most similar to.

This is research software, and its usefulness apart from its original context of detecting eukaryotes is not yet known :). Reference sequences are grouped by taxon, so its use with another reference database requires the provision of options `--refdb-format` or `--refdb-marker-to-taxon-path`. 


## Filtering options

Recommended presets are:

`" --min-read-mapq 30 --min-read-query-length 60 --min-read-match-identity 0.9 --min-taxon-num-markers 2"`
if using single best alignment per query.


If using multiple alignments, the following preset recommended if you're okay with relying on MCL clusters:
` --min-read-query-length 60 --min-taxon-num-markers 2 --min-taxon-num-reads 2 --min-taxon-better-marker-cluster-averages-ratio 1.01 --threshold-avg-match-identity-to-call-known-taxon 0.97  --threshold-num-taxa-to-call-unknown-taxon 1 --threshold-num-markers-to-call-unknown-taxon 4     --threshold-num-reads-to-call-unknown-taxon 8`

A simpler alternative is 
` --min-read-query-length 60 --min-taxon-num-markers 2 --min-taxon-num-reads 2 --min-taxon-fraction-primary-matches 0.5` 
but it does not deal with unknown taxa quite as well.

All filtering options are as follows:

| column | description |
| ------------- | ------------- | 
|`--min-read-mapq`                                   |when reading the input, skip alignments with MAPQ < min-read-mapq                                                                                                                                               |
|`--min-read-query-length`                           |when reading the input, skip alignments shorter than min-read-query-length                                                                                                                                      |
|`--min-read-match-identity`                         |when reading the input, skip alignments where the proportion of matching bases in the alignment is less than min-read-match-identity                                                                            |
|`--min-taxon-num-markers`                           |Only keep taxa with at least min-taxon-num-markers markers                                                                                                                                                      |
|`--min-taxon-num-reads`                             |Only keep taxa with at least min-taxon-num-reads reads                                                                                                                                                          |
|`--min-taxon-num-alignments`                        |Only keep taxa with at least min-taxon-num-alignments alignments                                                                                                                                                          |
|`--min-taxon-fraction-primary-matches`              |Only keep taxa where no more than min-taxon-fraction-primary-matches fraction of alignments is inferior / secondary                                                                                             |
|`--min-taxon-better-marker-cluster-averages-ratio`  |Only keep taxa where the ratio between markers which have at least average match identity relative to their clusters and markers with identity below average is at least min-taxon-better-cluster-averages-ratio|
|`--threshold-avg-match-identity-to-call-known-taxon`|Threshold on average match identity to return taxon in reference                                                                                                                                                |
|`--threshold-num-reads-to-call-unknown-taxon`       |To positively identify an unknown taxon (fits all criteria except match identity) expect this many reads from a taxon cluster                                                                                   |
|`--threshold-num-markers-to-call-unknown-taxon`     |To positively identify an unknown taxon (fits all criteria except match identity) expect this many markers from a taxon cluster                                                                                 |
|`--threshold-num-taxa-to-call-unknown-taxon`     |To positively identify an unknown taxon (fits all criteria except match identity) expect this many taxa from a taxon cluster                                                                                 |
### Reasons to apply filters

1. Very short alignments do not convey useful information
Our ancestors had to make do with 35-40bp shotgun reads, but we have longer ones - game changer for metagenomics! Still, a 100bp read can match on the last twenty bases at the end of a reference sequence (clipped alignments) or you could have configured the aligner to do local alignments instead of end-to-end. Either way, `--min-read-query-length` being something high enough (60 from EukDetect seems to work fine) addresses this problem.

2. Low identity matches are not taxon specific
An unknown species will match as a mixture of results. The clustering option `--min-taxon-better-marker-cluster-averages-ratio` tries to take care of removing the overall inferior evidence, and the `--threshold-avg-match-identity-to-call-known-taxon` only passes

The suggested value of 0.97 has been chosen empirically. Is a bit lower than CCMetagen's 0.9841 quoted from [Vu et al (2019)](https://pubmed.ncbi.nlm.nih.gov/29955203/), as this number was calculated from ribosomal subunits, we're not aware of a study that calculates average identity for BUSCOs. Most unknown taxa seem to match at around 0.9 identity, and a value 0.95 still permitted an unknown <i>Penicillinum</i> species to appear as a mixture.

3. Threshold of evidence for making claims
Claiming a eukaryote is present based on one read would be preposterous! It's not clear how many reads are "enough" to make a claim, and actually, no number of reads is enough because off-target matches follow patterns. We suggest gaining evidence from at least two markers, and a higher standard for ambiguous hits coming from species not in the reference. You can also only report unknown species if the results indicate its two nearest taxa with `--threshold-num-taxa-to-call-unknown-taxon` option.


### Other info

#### More output options
You can save an intermediate database produced by providing the `--sqlite-db-path` argument, and then query it with a `sqlite3` client.

#### Custom or different reference database
The default `--refdb-format` is `generic`, which tries to produce nice names, but may or may not match how you want it to. Set `--refdb-format` to `no-split` if you don't want the nice names, and if you want the taxa to be recognised really correctly, list a lookup table under `--refdb-marker-to-taxon-path`.

## Known issues
Quantitative information might be unreliable when there is very few reads.



For a large enough file, the sqlite query engine runs out of page numbers when doing a `group by`. In [my fork of HuMAnN with similar query code](https://github.com/wbazant/humann/commit/1dc767f855) I have solved this by adding `'PRAGMA max_page_count = 4294967292;'` before the `group by`. I've not yet ran into this issue when using this package.

## Credits
I took the method of splitting multiple aligned reads by a weighted average (with the second power of match identity as weights), and the method of calculating CPMs, from HuMAnN.
I was inspired by how MetaPhlAn calculates taxon CPMs from marker CPMs, although they have more options and I just ported the simple one.
I copied the package setup from EukDetect, and developed the package mostly in the context of alignments to the EukDetect reference.
An idea for what outputs might be useful to users comes jointly from these three tools.

For inspiration of what read properties are worth filtering on and how to do it, some credit goes to [TALON's `transcript_utils` file](https://github.com/mortazavilab/TALON/blob/master/src/talon/transcript_utils.py).

## How to cite

We now have a preprint on biorxiv: https://doi.org/10.1101/2022.03.09.483664 .

