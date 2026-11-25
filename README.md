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

## Installation

### Locally

You will need to [install Rust](https://www.rust-lang.org/tools/install) first.

Then you should be able to install combine-pbwt-sparsepaint by running:

```
cargo install --git https://github.com/mrcieu/combine-pbwt-sparsepaint
```

and the binary will be installed locally (in `$HOME/.cargo/bin` unless you have specified otherwise. Run `cargo help install` if you need to work out where):

```
~/.cargo/bin/combine -h
```

or clone the repository and build it:

```
git clone https://github.com/mrcieu/combine-pbwt-sparsepaint
cd combine-pbwt-sparsepaint
cargo build --release
```

and the binary gets built in the repo's directory:

```
❯ ls target/release
combine
```

### Docker

There is a Dockerfile in this repository which you can either use to make a distroless image containing `combine` (and nothing else), or an image with an alpine linux base (comment out / uncomment the relevant sections).

## Usage

### Input

Required:

1) `--chunkpathsfile`: a file containing the paths to the ouput files of `pbwt -paintSparse`, no header, one path per line.
2) `--snpcountsfile`: a file containing counts of the number of SNP present in each (vcf) file that was input to `pbwt -paintSparse`, no header, one count per line, in the same order as the outputs are given in `--chunkpathsfile`
3) `--maplengthsfile`: a file containing the map lengths for each file that was input to `pbwt -paintSparse`, no header, one value per line, in the same order as the outputs are given in `--chunkpathsfile`. Usually we run `pbwt` per-chromosome, so this is just a list of the map lengths for each chromosome according to the standard genetic map of your organism of choice.

Optional:

* `--restrict`: restrict the number of entries (lines) in the output to at most `2^31-1` when provided with no argument, by dropping the entries with the smallest values. This is the largest number of non-zero entries a sparse matrix can hold in R's [Matrix](https://cran.r-project.org/package=Matrix) package (and the maximum number of rows in a `data.table`), so the output can be loaded with e.g. `data.table::fread()` and `Matrix::sparseMatrix()`. Specific values can be provided: `--restrict=N` writes at most `N` entries. If several entries tie at the cut-off value, they are all dropped, so slightly fewer than `N` entries may be written. The equals sign `=` is required when using this flag with a value. Which entries are dropped is set by `--restrict-mode`.

* `--restrict-mode`: how `--restrict` drops entries, if there are more than `N`. Needs `--restrict`.
    * `global` (the default, and the only behaviour before v0.7.0): drop the entries with the smallest values in the whole matrix, using one cut-off for every row. Samples whose values are mostly small (e.g. those with fewer close matches in the data) can lose most or all of their entries.
    * `row-fraction`: keep the same fraction (`N` / the number of entries) of every row's entries: for a row with `n` entries, its largest `floor(n × N / total)`. Every sample keeps a share of its entries, so rows are affected evenly. As with `global`, ties at a row's cut-off are all dropped, and rounding down means slightly fewer than `N` entries are written.

* `--writesums`: write the sum of the scaled values in each row (i.e. for each `index1`) of the sparse matrix to `<PREFIX>.sums`. Sums are calculated before any entries are dropped by `--restrict`. If you use `--restrict`, this file will be written by default. See [Output](#output).

* `--threads`: number of threads to use for parallel processing of the output. Set to 8 by default. Must be `>=1`. The first part of the program (reading the data) is not parallelized, but the second part (writing the output) is parallelized, and speeds up with additional threads.

* `--out`: prefix for the output files. By default, this is set to `combined`, and the main output file will be called `<PREFIX>.txt.gz`. The sums file will be called `<PREFIX>.sums`, if written.

### Output

The main output is a gzip-compressed, space-separated text file with no header (`<PREFIX>.txt.gz`). It consists of three columns: `index1`, `index2`, and `value`, where `index1` is the (1-based) row index of the first sample in the original (vcf) file provided to `pbwt`, and `index2` is the (1-based) row index of the second sample. The `value` column contains the expected haplotype length, in map units, inherited between this pair of samples. E.g.:

```sh
❯ zcat combined.txt.gz | head -n4
1 11 10.50000
1 12 5.50000
1 540 11.50000
11 1 8.00000
```

The row sums file is a space-separated plain text file with no header, called `<PREFIX>.sums`. It has two columns, `index1` and `sum`. There is one line for every `index1` (sample) present in the input files, in increasing order, and `sum` is the total of that sample's values across all `index2`. E.g.:

```sh
❯ head -n2 combined.sums
1 7090.05071
2 7090.05072
```

The sums include every value `>= 0.000005` (the same values that are written to `<PREFIX>.txt.gz` without `--restrict`). With `--restrict`, they also include the values that were dropped, and samples whose values were all dropped are still listed. So you can use them to re-weight rows after filtering.

Each row sums to about `ploidy × total map length` (`ploidy` = 2 by default in `pbwt -paintSparse`), because `pbwt` adds both of a sample's haplotypes into one row. A sample whose sum is lower than that is missing from some of the input files (it had no matches to other samples on those chromosomes).

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

```sh
❯ combine \
	-c chunklength.files.txt \
	-s nsnps.txt \
	-m chr1-22.maplengths.txt \
	-r=100000 \
	--restrict-mode=row-fraction \
	-t 1 \
	-o chunklengths.restricted.rowfraction.100000
```

## Help

```
❯ combine -h
Usage: combine --chunkpathsfile <chunkpathsfile> --snpcountsfile <snpcountsfile> --maplengthsfile <maplengthsfile> [OPTIONS]

Options:
  -c, --chunkpathsfile <CHUNKPATHSFILE>
          A file containing the paths to the chunk files [default: ""]
  -s, --snpcountsfile <SNPCOUNTSFILE>
          A file containing the SNP counts for each chunk, in the same order as the chunk paths file [default: ""]
  -m, --maplengthsfile <MAPLENGTHSFILE>
          A file containing the map lengths for each chunk, in the same order as the chunk paths file [default: ""]
  -r, --restrict[=<RESTRICT>]
          The maximum number of entries (lines) to write to the output file. Which entries are dropped is set by --restrict-mode [Default with flag but no value: 2^31-1, the most R's Matrix package can hold]
      --restrict-mode <RESTRICT_MODE>
          How --restrict drops entries when there are too many: "global" drops the smallest values in the whole matrix; "row-fraction" keeps the same fraction of each row's entries (its largest values), so that rows with small values are not dropped disproportionately [default: global] [possible values: global, row-fraction]
  -w, --writesums
          Write the sum of each row (index1) to <OUT>.sums, as "index1 sum", counting all values >= 0.000005 including any dropped by --restrict. Always written if --restrict is used
  -t, --threads <THREADS>
          The number of threads to use for writing [default: 8]
  -o, --out <OUT>
          The prefix for the output file(s) [default: combined]
  -l
          Print a summary of licence information
      --licences
          Print full licence information, including all third-party notices
  -v, --version
          Print version
  -h, --help
          Print help (see more with '--help')
```
