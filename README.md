# Newick Format File Parser

Parsing of Newick tree format files written in R. Constructs a list of edges, leaf labels, and internal node labels and return the order in which leafs need to be aligned (progressive alignment).

## Dependencies (libraries)

* `stringr` : string operations

## Usage

```Rscript --vanilla newick_parse.R file.newick```
