# marker alignments
This a Python package to process and summarise alignments of metagenomic sequencing reads to reference databases of marker genes.

Alignments to marker genes are frequently used in shotgun metagenomics to quantify the abundance of taxa within a sample. Different methods of analysing the results are used in different software tools.

`marker_alignments` is a versatile and memory-efficient tool for working with alignments to marker genes. It uses a package `pysam` to read the files, so all common alignment formats should be supported. In addition to a few different outputs, you can also keep an intermediate database, and query it however you like.

## Install
```
pip3 install marker_alignments
```

## Usage
This is the most basic program:

```
process_marker_alignments --input tests/data/example.sam --output /dev/stdout
```

It produces a coverage report for each reference in the alignments file.

With information about reference format and number of reads used for the alignment, we can additionally attempt to quantify different taxa using the markers.

Here is how we can align a FASTA to EukProt, a reference database bundled with [EukDetect](https://github.com/allind/EukDetect), then summarise it at taxon level:

```

bowtie2 --omit-sec-seq --no-discordant --no-unal \
  -x ncbi_eukprot_met_arch_markers.fna \
  -1 ERR2749179_1.fastq \
  -2 ERR2749179_2.fastq \
  -S ERR2749179-eukprot.sam 

process_marker_alignments --input ERR2749179-eukprot.sam --output ERR2749179-eukprot-taxa.tsv \
  --refdb-format eukprot \
  --output-type taxon_all \
  --num-reads $(grep -c '^@' ERR2749179_1.fastq) \
```

For a larger alignment file, perhaps coming from a reference database such as ChocoPhlAn( a reference database bundled with [MetaPhlAn](https://github.com/biobakery/MetaPhlAn/) be sure to provide the `--sqlite-db-path` option with somewhere with plenty of disk space - then the program won't need to keep the alignments memory. 

You can query the database saved in this file with a `sqlite` client.

You can also integrate a package into a Python program:
```
import pysam
import re
from marker_alignments.main import read_alignments

alignment_store = read_alignments(
  alignment_file = pysam.AlignmentFile("ERR2749179.sam"),
  sqlite_db_path = ":memory:",
  pattern_taxon = re.compile("(^)"),
  pattern_marker = re.compile("(.*)")
)
markers = [ marker for marker in alignment_store.query('select distinct marker from alignment order by -coverage')]
```
### Custom refdb
The default `--refdb-format` is `generic`, which may or may not match how you want it to.
## Known issues
Quantitative information obtained by algining to EukProt might be a little bit dodgy since there are typically very few reads, and EukDetect was not designed for that.


For a large enough file, the sqlite query engine runs out of page numbers when doing a `group by`. In [my fork of HuMAnN with similar query code](https://github.com/wbazant/humann/commit/1dc767f855) I have solved this by adding `'PRAGMA max_page_count = 4294967292;'` before the `group by`. I've not yet ran into this issue when using this package.

## Credits
I took the method of splitting multiple aligned reads by a weighted average (with the second power of match identity as weights), and the method of calculating CPMs, from HuMAnN.
I was inspired by how MetaPhlAn calculates taxon CPMs from marker CPMs, although they have more options and I just ported the simple one.
I copied the package setup from EukDetect.
An idea for what outputs might be useful to users comes jointly from these three tools.
