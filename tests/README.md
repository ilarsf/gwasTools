# Tests

This repo uses lightweight smoke tests rather than a package-level test suite.

Run:

```sh
scripts/smoke_tests.sh
```

The script uses synthetic fixtures from `tests/fixtures/`, runs each `Rscript --help` command, generates QQ, Manhattan, binary power, and quantitative power plots in a temporary directory, and checks that each expected PNG exists and is non-empty.

It also checks adversarial synthetic cases:

- p-values outside `0 < P <= 1` fail clearly.
- declared MAF columns reject values greater than 0.5.
- Manhattan defaults accept the `CHROM` column used by fixtures.
- legacy `CHR` columns still work when `CHROM` is absent.
- mixed chromosome encodings such as `1`, `chr1`, `X`, and `23` are normalized.
- unsupported contig labels fail clearly.
- nearest-gene labeling reads the ignored `data/gene_positions/` cache by default and does not create or download references unless explicitly enabled.
- invalid power-plot sample-size vectors fail clearly.
- invalid quantitative power degrees of freedom fail clearly.

The fixtures are intentionally tiny and non-sensitive. Do not add real GWAS, genetic, phenotype, clinical, or PHI-like data to this directory.
