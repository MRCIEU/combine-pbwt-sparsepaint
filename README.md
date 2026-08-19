# combine

This is a Rust implementation of [combine_chunklength](https://github.com/YaolingYang/SparsePainter/tree/main/painting-pipeline/Compute%20haplotype%20components%20(HCs)), with some added functionality:

  1) accepts arguments from the command line
  2) restricts the output to `< 2^31` lines, necessary for loading it into R

## Usage

```
./combine -h
Usage: combine <chr_chunks_filenames_list> <chr_snp_counts_list> <chr_map_lengths_list> <nind> <outfile>
```
