# GEVA
Genealogical Estimation of Variant Age (GEVA)


This document describes the usage of the computational framework that we developed to estimate the age of the allele observed at a single locus in large-scale haplotype sample data.


> **Dating genomic variants and shared ancestry in population-scale sequencing data**  
> *Patrick K. Albers and Gil McVean*  
> doi: https://doi.org/10.1101/416610
###### See manuscript on bioRxiv: https://www.biorxiv.org/content/early/2018/09/13/416610

The GEVA method is described in detail in the Supplementary Text.



### Compilation ...
... is straightforward. Simply type `make` on the command line and the program compiles.
This generates a single executable called `geva_v1beta`.
You can use
```
./geva_v1beta --help
```
to see a list of available command line options.


## Conversion
In principle, the GEVA method operates on all haplotypes available in a given data set.
To avoid parsing the original source file every time GEVA is used, it parses the whole sample only once, which makes any subsequent loading of data into memory easier and much quicker.
The converstion creates three files; a binary file (`*.bin`), which contains the data, and two additional files (`*.marker.txt` and `*.sample.txt`), which list the parsed variant markers and samples, respectively.

Currently, only variant call format (VCF) files are supported, either uncompressed (`*.vcf`) or gzip compressed (`*.vcf.gz`).
Information about genetic distances is included during the conversion already; either by specifying a fixed recombination rate using the `--rec` option, or by providing a genetic map file using the `--map` option.
Recognized genetic map formats either have 3 or 4 columns, where the `Chromosome` column is optional; see example below.
```
Chromosome	Position(bp)	Rate(cM/Mb)	Map(cM)
chr1	55550	2.981822	0.000000
chr1	82571	2.082414	0.080572
chr1	88169	2.081358	0.092229
chr1	254996	3.354927	0.439456
...
```

The program operates on phased haplotype data, assuming diploid individuals.
It is required to convert data separately per chromosome; that is, to create a binary file for each chromosome.
Source files that combine data from multiple chromosomes need to be divided first, such that each input file contains data for one chromosome only.

To convert source file `DATA.vcf` (or `DATA.vcf.gz`) on the command line, use the examples provided below.
```
# fixed recombination rate, without a genetic map
./geva_v1beta --vcf DATA.vcf --rec 1e-8 --out NAME
```
or
```
# variable recombination rates, as provided through a genetic map
./geva_v1beta --vcf DATA.vcf --map /path/to/GENETIC_MAP_FILE --out NAME
```
The above creates the following files:
- `NAME.bin`
- `NAME.marker.txt`
- `NAME.sample.txt`

where `NAME` is the prefix specified using either the `-o` or `--out` argument.
Also, two additional files are created, a log file (`NAME.log`) and an error file (`NAME.err`). The latter is empty (0 bytes) if no errors or warnings were produced.
Note that `*.log` and `*.err` files are created in every run.


## Execution
The program loads the data contained in the generated `NAME.bin` file; specified using either the `-i` or `--input` argument.
The age estimation process relies on the detection of haplotype segments, shared between hundreds or throusands of haplotype pairs, which are detected relative to a given target position.
We developed a hidden Markov model (HMM) that uses empirically estimated emission and initial state probabilities, so as to be more robust towards data error.
The files that feed the HMM are provided in the `hmm` subdirectory, which are specified on the common line as given below.
```
--hmm ./hmm/hmm_initial_probs.txt ./hmm/hmm_emission_probs.txt
```
The program further requires the effective population size as a scaling parameter, using the `--Ne` argument, and the fixed rate of mutation (per site per generation), using the `--mut` argument.
Use the `--position` argument to specify the target variant (i.e. the chromosomal position of the variant whose age you want to estimate), or use the `--positions` argument to provide a batch file.
A batch file lists multiple target variants by position, where positions can be separated by any whitespace characters; e.g. newline, tab, or space.
See the examples given below.

```
# estimate allele age for the variant at position 12345678
./geva_v1beta -i NAME.bin -o RUN1 --position 12345678 --Ne 10000 --mut 1e-8 --hmm ./hmm/hmm_initial_probs.txt ./hmm/hmm_emission_probs.txt
```

```
# estimate allele age multiple variants, listed in file BATCH.txt
./geva_v1beta -i NAME.bin -o RUN1 --positions /path/to/BATCH.txt --Ne 10000 --mut 1e-8 --hmm ./hmm/hmm_initial_probs.txt ./hmm/hmm_emission_probs.txt
```

Note that the alternative allele is assumed to be the derived allele. The distribution of the derived allele in the sample is used to determine the pairing of haplotypes, i.e. to form concordant and discordant pairs.
Concordant pairs are pairs where both haplotypes carry the target allele, and discordant pairs consist of one carrier and one non-carrier.
The program samples up to a specified number of pairs from each group.
To change the default sampling limits (100 for both groups), you can use the `--maxConcordant` and `--maxDiscordant` command line options.
See example below.
```
./geva_v1beta -i NAME.bin -o RUN1 --positions /path/to/BATCH.txt --maxConcordant 500 --maxDiscordant 500 --Ne 10000 --mut 1e-8 --hmm ./hmm/hmm_initial_probs.txt ./hmm/hmm_emission_probs.txt
```

Again, use the `-o` or `--out` argument to specify the prefix for the files generated.


## Output
By executing the program as described above, two result files are created (plus a `*.log` and a `*.err` file):
- `RUN1.pairs.txt`
- `RUN1.sites.txt`

The **pairs** file contains the results of all pairwise analyses conducted to estimate the age of the variant(s) contained in the **sites** file.
The format of each file type is described below.

*(to be completed)*


## Comments

*(to be completed)*

