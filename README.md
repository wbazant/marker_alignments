# marker alignments
This a Python package to process and summarise alignments of metagenomic sequencing reads to a reference database of marker genes. You can use it in conjunction with an aligner like `bowtie2` to produce an estimate of taxa present in a metagenomic sample.

The package was developed in the context of looking for eukaryotes - most of the facilities are for producing good guesses from small amounts of potentially unreliable information. There are read level filters, clustering facilities for making sense of multiple alignments per query, and a number of thresholds.

## Installation
First clone this repository, then install the dependencies and this package:
```
git clone git@github.com:wbazant/marker_alignments.git
cd marker_alignments

pip3 install -r requirements.txt .
```

## Usage

### Basic example
```
marker_alignments --input tests/data/example.sam --output /dev/stdout
```

It produces a coverage report for each reference in the alignments file.


### Example (EukDetect reads)

Here is how we can align a FASTA to a reference database bundled with [EukDetect](https://github.com/allind/EukDetect), then summarise it at taxon level:

```

bowtie2 --omit-sec-seq --no-discordant --no-unal \
  -x ncbi_eukprot_met_arch_markers.fna \
  -1 ERR2749179_1.fastq \
  -2 ERR2749179_2.fastq \
  -S ERR2749179-eukprot.sam 

marker_alignments --input ERR2749179-eukprot.sam --output ERR2749179-eukprot-taxa.tsv \
  --refdb-format eukprot \
  --refdb-marker-to-taxon-path busco_taxid_link.txt \
  --output-type taxon_all \
  --num-reads $(grep -c '^@' ERR2749179_1.fastq) \
```

If this is your use case, try the Nextflow pipeline [wbazant/marker-alignments-nextflow](https://github.com/wbazant/marker-alignments-nextflow).


## How to use this package
The basic workflow is to give it an alignment file, and look at reports produced. Setting up various filters reduces noise enough that the taxonomic profile obtained can be passed on to other tools, but this is still research


### Filtering options
Recommended presets are:

`" --min-read-mapq 30 --min-read-query-length 60 --min-read-match-identity 0.9 --min-taxon-num-markers 2"`
if using single best alignment per query.


If using multiple alignments, the following options are recommended if you're okay with relying on MCL clusters:
` --min-read-query-length 60 --min-taxon-num-markers 2 --min-taxon-num-reads 4 --min-taxon-better-marker-cluster-averages-ratio 1.01 --threshold-avg-match-identity-to-call-known-taxon 0.97  --threshold-num-taxa-to-call-unknown-taxon 1 --threshold-num-markers-to-call-unknown-taxon 4     --threshold-num-reads-to-call-unknown-taxon 8`

A simpler alternative is 
` --min-read-query-length 60 --min-taxon-num-markers 2 --min-taxon-num-reads 4 --min-taxon-fraction-primary-matches 0.5` 
but it does not deal with unknown taxa quite as well.

All filtering options are as follows:
| ------------- | ------------- | 
|`--min-read-mapq`                                   |when reading the input, skip alignments with MAPQ < min-read-mapq                                                                                                                                               |
|`--min-read-query-length`                           |when reading the input, skip alignments shorter than min-read-query-length                                                                                                                                      |
|`--min-read-match-identity`                         |when reading the input, skip alignments where the proportion of matching bases in the alignment is less than min-read-match-identity                                                                            |
|`--min-taxon-num-markers`                           |Only keep taxa with at least min-taxon-num-markers markers                                                                                                                                                      |
|`--min-taxon-num-reads`                             |Only keep taxa with at least min-taxon-num-reads reads                                                                                                                                                          |
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
Claiming a eukaryote is present based on one read would be preposterous! It's not clear how many reads are "enough" to make a claim, and actually, no number of reads is enough because off-target matches follow patterns. We suggest EukDetect's rule of 4 reads in 2 markers, and doubling it for unknown species. You can also only report unknown species if the results indicate its two nearest taxa with `--threshold-num-taxa-to-call-unknown-taxon` option.


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

