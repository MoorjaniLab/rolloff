# ROLLOFF

ROLLOFF is a method for estimating the time of admixture. Details of the method and statistic can be found in Moorjani et al. 2013 
<a href="https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1001373">Reconstructing Roma History from Genome-Wide Data</a>. This method was first introduced in <a href="https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1001373">Moorjani et al. 2011</a> and distributed as part of <a href="https://github.com/DReichLab/AdmixTools"> ADMIXTOOLS </a>. The latest implementation is more reliable and robust to biases that can be introduced due to strong founder events that may postdate admixture.

#### Installation
We have placed source code for all C executables in the src/ directory, 
for users who wish to modify and recompile our programs.  For example, to
recompile the programs, type

```
cd rolloff
make clobber
make rolloff      
```

If you are building on a Mac, you will need gsl and openblas installed. 
```
brew install gsl
brew install homebrew/science/openblas
Uncomment the lines in src/Makefile that modify the CFLAGS and LDFLAGS. 
```

#### Input
ROLLOFF requires that the input data is available in one of these formats (See https://reich.hms.harvard.edu/software/InputFileFormats). To convert to the appropriate format, one can use CONVERTF program (See https://github.com/argriffing/eigensoft/tree/master/CONVERTF for details). 

#### Command line 
```
./rolloff -p $parfile >$logfile
```
$logfile: Name of the logfile. The ROLLOFF program prints various statistics to standard output which should be directed to the logfile.  <br />
$parfile: Name of parameter file.  <br />

#### Parameter file
```
genotypename:   input genotype file   # in eigenstrat format
snpname:   input snp file             # in eigenstrat format
indivname:   input indiv file         # in eigenstrat format
poplistname:   rolloff.ref.           # This file contains the names of the two ancestral populations. There is one population on each line.
admixlist:   rolloff.target.          # This file contains the name of the admixed population. There is one population on each line.
binsize:   binsize (in Morgans).      # Range is from 0-1. Optimal binsize of 0.001 is recommended.
maxdis:   maximum_distance (in Morgans). # Range is 0-1. For quicker runs, use max_distance < 1.0. However, for recent admixture, ensure that max_distance is greater than the expected admixture LD blocks.
seed:   rand_num.                     #Random seed to ensure reproducibility of runs. 
runmode:   1                          # internal parameter
chithresh:                            # chi-square threshold (default: 6.0). 
checkmap: YES/NO.                     # Checks if physical positions are correlated to genetic position. Error, if very high correction. 
```

##### Optional paramaters
```
weightname:   weight_file.        # Contains a weight for each SNP to be included in the run. If this parameter is not specified, the program uses the allele frequency differentiation between the ancestral populations as the weight for each SNP. 
mincount:   number.               # Contains the minimum number of admixed individuals required for the run. Default = 5.
minparentcount: number.           # Contains the minimum number of ancestral individuals from each population that is required for the run. Default = 10.
chrom: chromosome_number.         # The analysis is limited to the specified chromosome only.
nochrom:   chromosome_number.     # The specified chromosome is excluded from the analysis.
admixlist: admix_list.            # If you want run ROLLOFF for a list of populations, use flag admixlist to specify a list of populations. There is one population on each line. NOTE: The ancestral populations remain the same.
badsnpname:    badsnp_list.       # File contains a list of SNPs to be excluded from the analysis. 
ransample:   random_number.       # Program will run rolloff with a random set (n = random_number) of samples from the admixed population.
numchrom:	integer_number       # the count of chromosomes of your organism
```

#### Output
The program generates several output files (and some tempfiles).
```
output.out                        # This file contains output for the entire genome. The estimates of covariance values at various genetic distances, binned according to input values.
output.out:$chr, where chr=1-22   # These files contain the output for the jackknife where we remove one chromosome ($chr) in each run. 
```

Also read README.ROLLOFF_OUTPUT for fitting an exponential to the output of Rolloff and inferring the standard error using jackknife.

#### Example 
An example run is available in ``src/example/`` directory. 

#### Support
For questions, please email Priya Moorjani (moorjani@berkeley.edu)

