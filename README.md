# Set of scripts to summarize vcf results.

## INSTALLATION

Clone the repo:


```bash

git clone git@github.com:sschmeier/vcfcompile.git

```

### Requirements


 - Python 3
 - Otherwise nothing special. Uses only standard libs for now.



## vcfcompile

### DESCRIPTION

Simple script to read a bunch of vcf-files with SNPs and find the list of unique SNPs.
For each variant, extract for each file the variant's QD or QUAL value and put in table form.
Gives an overview of overlapping variants and quality values in different samples.

Prints to standard out. Some stats go to standard error.

Untested on very large vcf-files. In the future need to implement cyvcf for speed. Right now its used for filtered SNPs.

### Usage

```bash
python vcfcompile.py --snpeff data/*.vcf(.gz) > table.txt
```

### Output

| CHROM | POS      | ID        | REF | ALT | GENES            | FILE1.vcf.gz | FILE2.vcf.gz | ... |
|-------|----------|-----------|-----|-----|------------------|--------------|--------------|-----|
| chr17 | 16382069 | rs1060079 | T   | C   | UBB:HIGH;UBB:LOW | 2.99         | 3.64         | ... |
| ...   |          |           |     |     |                  |              |              |     |



## vcfSetStats.py

### DESCRIPTION

For a vcf-file that was compiled with `gatk3 CombineVariants`. 
The idea here is that the same sample was processed by different callers and the vcf is the combined file, we can use this script to investigate numbers number of variants called by any combination of callers.

### Usage

```bash
python vcfSetStats.py file.vcf.gz > table.tsv
```

### Output

A table with caller combination, number of callers, number of variants called, pct of variants called.



## TODO

 - Make use of cyvcf (https://github.com/arq5x/cyvcf) for speed.



## LICENCE

MIT, 2018-2019, copyright Sebastian Schmeier
s.schmeier@gmail.com // https://www.sschmeier.com
