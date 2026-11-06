# Command-line entry point for batch locus analysis

The input is an RDS list containing the named arguments of
[`susieout_fit()`](https://weinstocklab.github.io/SuSiEOut/reference/susieout_fit.md).
The output is a complete fit RDS, with an optional tab-separated support
file. The installed script `exec/susieout` calls this function.

## Usage

``` r
susieout_cli(args = commandArgs(trailingOnly = TRUE))
```

## Arguments

- args:

  Character vector of command-line arguments: `--input`, `--output`, and
  optionally `--summary`, `--L`, `--lambda-total`, and `--out-fraction`.
  Use `--help` for usage.

## Value

Invisibly, `0` on success. The fit and optional support table are
written to the paths given.

## Examples

``` r
susieout_cli("--help")
#> Usage: susieout --input LOCUS.rds --output FIT.rds [--summary SUPPORT.tsv] [--L 10] [--lambda-total 1] [--out-fraction 0.10]
```
