# combine

Combine the (usually) per-chromosome output of [pbwt](https://github.com/richarddurbin/pbwt) `-paintSparse` into a single file, re-scaling chunk lengths to be in map units rather than number of SNPs.

For use when calculating [haplotype components](https://www.nature.com/articles/s41467-025-57601-3).

The software in this repository is a re-implementation of [combine_chunklength](https://github.com/YaolingYang/SparsePainter/tree/main/painting-pipeline/Compute%20haplotype%20components%20(HCs)).

### Haplotype Component Calculation Pipeline

Documentation for the original haplotype component calculation pipeline [can be found here](https://github.com/YaolingYang/SparsePainter/tree/main/painting-pipeline/Compute%20haplotype%20components%20(HCs)). Briefly, the steps are:

1) (phased) VCF -> per-chromosome sparse matrices of haplotype sharing lengths (in number of SNPs) using [pbwt](https://github.com/richarddurbin/pbwt)
2) per-chromosome sparse matrices of haplotype sharing lengths (in number of SNPs) -> total sparse matrix of haplotype sharing lengths (in map units)
3) sparse matrix of haplotype sharing lengths (in map units) -> haplotype components (HCs) using singular value decomposition

The program in this repository performs step 2.

## Usage

### Input

Required:

1) `--chunkpathsfile`: a file containing the paths to the ouput files of `pbwt -paintSparse`, no header, one path per line.
2) `--snpcountsfile`: a file containing counts of the number of SNP present in each (vcf) file that was input to `pbwt -paintSparse`, no header, one count per line, in the same order as the outputs are given in `--chunkpathsfile`
3) `--maplengthsfile`: a file containing the map lengths for each file that was input to `pbwt -paintSparse`, no header, one value per line, in the same order as the outputs are given in `--chunkpathsfile`. Usually we run `pbwt` per-chromosome, so this is just a list of the map lengths for each chromosome according to the standard genetic map of your organism of choice.

Optional:

* `--restrictrows`: restrict the number of rows in the output (to `< 2^31` when provided with no argument). Specific values can be provided: `--restrictrows=N` restricts the number of rows to `N`. The equals sign `=` is required when using this flag with a value.

* `--writerowsums`: write the sum of the scaled values for each row in the sparse matrix (prior to any filtering because of `--restrictrows`) to a file. If you use `--restrictrows`, this file will be written by default.

* `--threads`: number of threads to use for parallel processing of the output. Set to 8 by default. Must be `>=1`. The first part of the program (reading the data) is not parallelized, but the second part (writing the output) is parallelized, and speeds up with additional threads.

* `--out`: prefix for the output files. By default, this is set to `combined`, and the main output file will be called `<PREFIX>.txt.gz`. The rowsums file will be called `<PREFIX>.rowsums`, if written.

### Output

The main output is a gzip-compressed tab-separated text file with no header (`<PREFIX>.txt.gz`). It consists of three columns: `index1`, `index2`, and `value`, where `index1` is the (1-based) row index of the first sample in the original (vcf) file provided to `pbwt`, and `index2` is the (1-based) row index of the second sample. The `value` column contains the expected haplotype length, in map units, inherited between this pair of samples. E.g.:

```sh
❯ zcat combined.txt.gz | head -n3
1     11      10.5
1     12      5.5
1     540     11.5
11    1       8.0
```

The rowsums file is a plain text file called `<PREFIX>.rowsums`. It contains the sum of the values in the value column for each row index. E.g.:

```sh
❯ head -n1 combined.rowsums 
27.5
```

### Example command line use

```sh
❯ combine \
	-c chunklength.files.txt \
	-s nsnps.txt \
	-m chr1-22.maplengths.txt \
	-o chunklengths
```

```sh
❯ combine \
	-c chunklength.files.txt \
	-s nsnps.txt \
	-m chr1-22.maplengths.txt \
	-r \
	-t 1 \
	-o chunklengths.restricted
```

```sh
❯ combine \
	-c chunklength.files.txt \
	-s nsnps.txt \
	-m chr1-22.maplengths.txt \
	-r=100000 \
	-t 1 \
	-o chunklengths.restricted.100000
```

## Help

```
❯ combine -h
Usage: combine [OPTIONS] --chunkpathsfile <CHUNKPATHSFILE> --snpcountsfile <SNPCOUNTSFILE> --maplengthsfile <MAPLENGTHSFILE> --nsample <NSAMPLE>

Options:
  -c, --chunkpathsfile <CHUNKPATHSFILE>
          A file containing the paths to the chunk files
  -s, --snpcountsfile <SNPCOUNTSFILE>
          A file containing the SNP counts for each chunk
  -m, --maplengthsfile <MAPLENGTHSFILE>
          A file containing the map lengths for each chunk
  -r, --restrictrows[=<RESTRICTROWS>]
          The maximum number of rows to write to the output file [Default with flag but no value: 2^31-1]
  -w, --writerowsums
          Write the row sums (prior to any dynamic filtering). If --restrictrows is in effect, this file will be written anyway
  -t, --threads <THREADS>
          The number of threads to use for writing [default: 8]
  -o, --out <OUT>
          The prefix for the output file(s) [default: combined]
  -h, --help
          Print help
  -V, --version
          Print version
```
