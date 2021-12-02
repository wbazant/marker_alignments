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

| Parameter  | Description |
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


