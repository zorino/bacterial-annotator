# Bacterial-Annotator GEM

Ruby GEM to annotate bacterial genomes based on a reference genome and also complete the annotation with external or remote databases.

To install the gem locally with this package
```shell
gem build bacterial-annotator.gemspec
gem install bacterial-annotator-0.0.1.gem

bacterial-annotator -h
```

Or you could simply install it from rubygems.org with that command
```shell
gem install bacterial-annotator
```

The blat aligner and the prodigal ORF finder will be installed as dependencies if you accept their licenses.

## Example of usage


* Bacterial genome annotation
``` shell
bacterial-annotator annotate -i my_genome.fasta -g reference_genome.gbk -o annotation_output_dir --pidentity 80 --pcoverage 90

# reference_genome.gbk could also only be the genome ID which will be fetched from NCBI

```

* Bacterial genome comparison
``` shell
# in a folder where you want to compare diffente genome annotations
bacterial-annotator compare --refgenome reference_genome.gbk --pidentity 80 --align prot --phylogeny --software fasttree *.gbk
```
