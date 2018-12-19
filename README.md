# vcfcompile

## DESCRIPTION

Simple script to read a bunch of vcf-files with SNPs and find the list of unique SNPs.
For each variant, extract for each file the variant's QD or QUAL value and put in table form.
Gives an overview of overlapping variants and quality values in different samples.

Prints to standard out. Some stats go to standard error.

Untested on very large vcf-files. In the future need to implement cyvcf for speed. Right now its used for filtered SNPs.


## INSTALLATION

### Requirements

 - Python 3
 - Otherwise nothing special. Uses only standard libs for now.


## USAGE

```bash
python vcfcompile.py data/*.vcf(.gz) > table.txt
```


## TODO

 - Make use of cyvcf (https://github.com/arq5x/cyvcf) for speed.


## VERSION HISTORY

- 0.0.1    2018    Initial version.


## LICENCE

MIT, 2018, copyright Sebastian Schmeier
s.schmeier@gmail.com // https://www.sschmeier.com
