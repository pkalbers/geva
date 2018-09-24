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

### pairs
The `*.pairs.txt` file has the following fields:
- `MarkerID` : Focal variant; internal ID of the marker (as given in `*.marker.txt` file).
- `Clock` : The clock model used; mutation clock (`M`), recombination clock (`R`), joint clock (`J`).
- `SampleID0` : Internal ID of the first individual (diploid!) in the pair (as given in `*.sample.txt` file).
- `Chr0` : Indicator of the haplotype in the first individual; either `0` or `1` (first or second haplotype).
- `SampleID1` : Same as `SampleID0` but for the second individual in the pair.
- `Chr1` : Same as `Chr0` but for the second individual in the pair.
- `Shared` : Is pair concordant (`1`) or discordant (`0`) ?
- `Pass` : Did this pair pass quality control, using a heuristic method for rejecting pairs?
- `SegmentLHS` : Breakpoint detected (HMM) on *left* hand side from the focal site; given as internal ID of marker.
- `SegmentRHS` : Breakpoint detected (HMM) on *right* hand side from the focal site; given as internal ID of marker.
- `Shape` : Value of *shape* paramter; Gamma distirbution.
- `Rate` : Value of *rate* paramter; Gamma distirbution.

The *shape* and *rate* parameters are used to obtain a posterior distribution on the TMRCA of a given pair. Allele age is estimated from the composite posterior distribution, which combines the pairwise TMRCA posteriors available for a given focal variant.

### sites
The `*.sites.txt` file has the following fields:
- `MarkerID` : Focal variant; internal ID of the marker (as given in `*.marker.txt` file).
- `Clock` : The clock model used; mutation clock (`M`), recombination clock (`R`), joint clock (`J`).
- `Filtered` : Was allele age computed before (`0`) or after (`1`) quality control (heuristic filtering of pairs)?
- `N_Concordant` : The number of available concordant pairs (before or after filtering).
- `N_Discordant` : The number of available discordant pairs (before or after filtering).
- `PostMean` : The *mean* of the composite posterior distribution.
- `PostMode` : The *mode* of the composite posterior distribution.
- `PostMedian` : The *median* of the composite posterior distribution.

Note that all coordinates refer to the internally used IDs that are given in the `*.marker.txt` and `*.sample.txt` files generated during compilation.  

Several results are given for each focal variant; there is one allele age estimate for each clock model, first, based on all pairs analysed and, second, based on the set of pairs retained after quality control.
This is distinguished by the `Filtered` field in the `*.sites.txt` file, and by the `Pass` field in the `*.pairs.txt` file.

However, the exact heuristic filtering algorithm used here differs from the one that we used for generating the results presented in the paper. We applied an external script to filter pairs in the `*.pairs.txt` file, from which we estimated allele age.  The script runs in R and is provided: `estimate.R`  
To run this script, for example on the generated `RUN1.pairs.txt` file, execute on the command line
```
Rscript estimate.R /path/to/RUN1.pairs.txt 10000
```
where `10000` refers to the scaling parameter, Ne.  
The above creates a new "sites" file, but now named `RUN1.sites2.txt`.


## Comments
The GEVA framework, as it is currently implemented, has a few known bugs; listed below.

- Estimating the age of the very first or the very last variant of a chromosome may produce `Segmentation fault` errors, sometimes.
- The formation of concordant pairs may sometimes fail if there is no heterozygous state found in the sample for a given target variant.
- Memory allocation grows exponentially with allele frequency.

Due to the latter, it is highly recommended to keep batch files small, in the order of hundreds.

