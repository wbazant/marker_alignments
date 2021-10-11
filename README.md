# marker alignments
This a Python package to process and summarise alignments of metagenomic sequencing reads to a reference database of marker genes. You can use it in conjunction with an aligner like `bowtie2` to produce an estimate of taxa present in a metagenomic sample.


## Installation
First clone this repository, then install the dependencies and this package:
```
git clone git@github.com:wbazant/marker_alignments.git
cd marker_alignments

pip3 install -r requirements.txt .
```

## Usage
This is the most basic program:

```
marker_alignments --input tests/data/example.sam --output /dev/stdout
```

It produces a coverage report for each reference in the alignments file.

With information about reference format and number of reads used for the alignment, we can additionally attempt to quantify different taxa using the markers.

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

For a larger alignment file, perhaps coming from a reference database such as ChocoPhlAn( a reference database bundled with [MetaPhlAn](https://github.com/biobakery/MetaPhlAn/) be sure to provide the `--sqlite-db-path` option with somewhere with plenty of disk space - then the program won't need to keep the alignments memory. 

You can query the database saved in this file with a `sqlite` client.

You can also integrate a package into a Python program:
```
import pysam
import re
from marker_alignments.summarize.main import read_alignments

alignment_store = read_alignments(
  alignment_file = pysam.AlignmentFile("ERR2749179.sam"),
  sqlite_db_path = ":memory:",
  pattern_taxon = re.compile("(^)"),
  pattern_marker = re.compile("(.*)")
)
markers = [ marker for marker in alignment_store.query('select distinct marker from alignment order by -coverage')]
```

## Filtering the results


| Parameter  | Filtering level | Description |
| ------------- | ------------- | ------------ |
| --min-read-mapq  | aligned read   | MAPQ field in the alignment file |
| --min-read-query-length  | aligned read | Length of the alignment | 
| --min-read-match-identity  | aligned read  | Fraction of the alignment with matching bases |
| --min-taxon-num-markers | taxon | Number of markers where the taxon had aligned reads |
| --min-taxon-num-reads | taxon | Number of reads aligned to taxon's markers |

### Custom refdb
The default `--refdb-format` is `generic`, which tries to produce nice names, but may or may not match how you want it to. Set `--refdb-format` to `no-split` if you don't want the nice names, and if you want the taxa to be recognised really correctly, list a lookup table under `--refdb-marker-to-taxon-path`.

## Known issues
Quantitative information obtained by aligning to a EukDetect reference database might be a little bit dodgy since there are typically very few reads.

The results might contain a lot of noise, and excluding low-length and low quality alignments might help. Passing the alignments file through `samtools view -m40 -q30 -h` before summarizing might go a long way.


For a large enough file, the sqlite query engine runs out of page numbers when doing a `group by`. In [my fork of HuMAnN with similar query code](https://github.com/wbazant/humann/commit/1dc767f855) I have solved this by adding `'PRAGMA max_page_count = 4294967292;'` before the `group by`. I've not yet ran into this issue when using this package.

## Credits
I took the method of splitting multiple aligned reads by a weighted average (with the second power of match identity as weights), and the method of calculating CPMs, from HuMAnN.
I was inspired by how MetaPhlAn calculates taxon CPMs from marker CPMs, although they have more options and I just ported the simple one.
I copied the package setup from EukDetect, and developed the package mostly in the context of alignments to the EukDetect reference.
An idea for what outputs might be useful to users comes jointly from these three tools.

For inspiration of what read properties are worth filtering on and how to do it, some credit goes to [TALON's `transcript_utils` file](https://github.com/mortazavilab/TALON/blob/master/src/talon/transcript_utils.py).


## Explanation of the method

### Definitions
**taxon** - an organism, in this context an organism that had its genome sequenced, appears in the reference database, and may or may not be present in the sequenced sample

**marker** - a DNA sequence of a gene in a taxon that is assumed to be unique to the taxon

**BUSCO** - a family of genes that are mostly present in each taxon and mostly single-copy

**reference database** - one of the inputs for an aligner, in this context it's a reference of markers that can be matched to

**alignment / hit** - a read in the sequenced sample that was found by an aligner to match a marker in the reference

Aligners like bowtie2 can report multiple alignments per read. In that case, we can distinguish:
**primary alignment** - the best match for the read (based on alignment score, sequence identity, or other metric)

**secondary alignments** - alignments corresponding to matches for a read that are not the best one


We would also like to talk about kinds of hits:

**true hit** - a read that matches a marker M for a taxon A because of A being present in the sample 

**false hit** - a hit that is not a true hit


### Sources of false hits
Hits could happen for no biological reason - the sequencing process introduces errors, short sequences can coincide by chance, and there are common genetic sequences that are under high selective pressure to not vary, like common binding sites. Using good enough reads and long enough matches should minimise these kinds of hits.

Additionally, the process of producing alignments could have an error rate as well, because aligners like `bowtie2` use approximate methods. However, they work well - we can assume they report the best match possible for each read, and when a match exists they will report it.

While false hits are inevitable, true hits are, actually, impossible in metagenomics. It is unlikely that an organism present in a sample will be exactly like one of the taxa in the reference, unless the sample is exactly what was used to create the sequence. Instead, we will allow ourselves a notion of an evolutionary past - organisms being related through their common ancestors. This lets us produce a following dichotomy:

**target hit** - a read that comes from an orgamism A, and matches a marker M' for a taxon A', such that A' is the closest taxon to A

**off-target hit** - a read that comes from an organism A, but matches a marker M' for a taxon A'', where A, A' and A'' are related species but A is closer to A' than A''

#### A model for off-target hits
Suppose an organism A has a version b(A) of a BUSCO b, and the reference contains markers b(A1)... b(An) for taxa A1... An. Let us say A is most similar to A1perhaps it's another strain of the same species. Assume also a least common ancestor A0 of A and A1, and A00 of A0, A2,...An.
As mutations accumulate over time, we can predict b(A) will be most similar to b(A1), but - in places where A1 has diverged from A0 - some segments of b(A) are most similar to other b(Ai). Some segments of b(A) could also be equally similar in all b(Ai), if there has been reason for that sequence not to change since the joint common ancestor A00.

This generates a prediction of what to expect when reads from b(A) are aligned to each of the b(Ai) - match identity should form a distribution, with increased evolutionary distance between A and Ai leading to lower average identity and fewer matches. Additionally, competitive alignment of reads from b(A) between b(A1)...b(An) will not entirely favour the closest A1.

If we ask an aligner to report a single best alignment for each read, we expect to see h1 hits to b(A1), and smaller amounts hi of hits for other b(Ai), such that H = sum(h1...hn) is proportional to the count of reads coming from sequencing b(A). We also expect ratios of H to hi to be related to sequence similarity of b(A) and b(Ai) - the further A is from A1, the larger the number of off-target hits H - h1.

The effect of b(A1)...b(An) 'competing' for the best alignment of each read is illuminated when an aligner is asked to report all reads. Some of h2..hn are then accompanied by secondary alignments to b(A1) - call them s1, and some of h1 will be accompanied by secondary alignments to b(A2)...b(An), s2...sn. If h1 is much larger than hi, s1 should be much smaller than si and thus h1/s1 should be larger than hi/si and independent of H.

Thus secondary alignments help us differentiate the presence of an organism A reported as many hits to A1 and fewer hits to A2 from the presence of two unrelated organisms X and Y. For example, we can report a ratio of primary to secondary alignments for each taxon.

It is possible for b(A1) to be very different from b(A), or missing from the reference entirely. Perhaps the genome of A1 is incorrectly annotated, or A1 has lost b when adapting to its niche.


#### Example 1 - off-target hits in more than just a genus
SRR6262267 is a run from a sample dominated by *Trichosporon asahii* - according to SRA Traces, 23.68% of the reads in the sample can be attributed to this organism ([source](https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR6262267)).

TODO an example - perhaps a heatmap of counts for BUSCOs in each Trichosporon?
Actually, EukDetect reports only T. asahii, because the other off-target hits are for the same genus - find something :).
TODO this will require some visualisation tools.

#### Example 2 - sticking to the reference too closely brings up nothing
Mucor example, demonstrate EukDetect returns nothing which it really should.


### Marker clusters
- empirically, using common alignments
- how we build a graph and run MCL on it
- how the competitive mapping is realised through "winning" 



