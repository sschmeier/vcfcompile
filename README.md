# vcfcompile

## DESCRIPTION

Simple script to read a bunch of vcf-files and compile a table of all variants and extract for each file the QD value of the SNPs and list them for each file. Prints to standard out. Some stats go to standard error.

Untested for very large vcf-files. IN the future need to implement cyvcf for speed.


## INSTALLATION

Nothing special. Uses only standard libs.


## USAGE

```bash
python vcfcompile.py *.vcf(.gz) > table.txt
```


## TODO

 - Make use of cyvcf for speed.


## VERSION HISTORY

- 0.0.1    2018    Initial version.


## LICENCE

MIT, 2018, copyright Sebastian Schmeier
s.schmeier@gmail.com // https://www.sschmeier.com
