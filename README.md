# go50k – Barley 50K SNP Processing

**go50k** (*Am Fear Liath Mòr*) is a command-line pipeline for processing Illumina
Infinium Barley 50K SNP array data. It takes the raw GenomeStudio *Final Report*
exports, reshapes them, clusters the genotypes with a parallelised R script and
converts the resulting AB calls into nucleotide (ACGT) calls.

Developed by Paul Shaw (ICS), Jim McNicol (BioSS) and Malcolm Macaulay (CMS),
International Barley Hub / James Hutton Institute, Invergowrie, Scotland, DD2 5DA.

---

## What it does

The program runs a series of stages, printing progress to the terminal:

| Stage | Description | Code |
|-------|-------------|------|
| 1/6 | Reads `Project.csv` and builds a map of Sentrix ID (`column 7_column 8`) to `Sample_ID` | [projectfile](projectfile/loadprojectfile.go) |
| 2/6 | Reads `submission_template.txt` | [submissionfile](submissionfile) |
| 4/6, 5/6 | Reads `FinalReportR.txt` (intensities) and `FinalReportTheta.txt` (theta values) **in parallel** using goroutines and transposes each into a marker × sample CSV named `OUT_<original name>` | [loaddata](loaddata/loaddata.go) |
| 6/6 | Runs [cluster50kparallel.R](cluster50kparallel.R) via `Rscript` to call genotypes (k‑means clustering) | [main.go](main.go) |
| post | Converts AB calls to ACGT using the `calls.txt` reference, then normalises odd calls and writes allele statistics | [main.go](main.go) |

### Clustering (R)

[cluster50kparallel.R](cluster50kparallel.R) is a parallelised, `data.table`-based
version of Jim McNicol's k‑means clustering method. For every SNP it:

1. Removes samples with missing theta/intensity values.
2. Identifies *null* (low intensity) genotypes.
3. Runs k‑means with 1–5 clusters on theta and chooses the number of clusters
   (1–3) using thresholds on the reduction in within‑cluster sum of squares and a
   minimum distance between cluster means.
4. Labels clusters `AA`, `AB`, `BB` (or ` N` for nulls) and records cluster
   statistics (sizes, means, SDs, separation, compactness, flag).

SNPs are distributed across all CPU cores but one using `foreach`/`doParallel`.
The clustering parameters (`my.cut.12`, `my.mindist`, `my.seed`, …) are defined
near the bottom of the R script.

### Post-processing (Go)

* **AB → ACGT** (`convertAB2ACGT`): reads `sample_data/ids.csv`, looks up each
  marker's A/B alleles in `calls.txt` and writes `sample_data/final_calls.txt`.
  Calls that cannot be translated are written as `-`/`--`.
* **Normalisation** (`ProcessGenotypeData`): cleans up ambiguous values
  (e.g. `-A` → `AA`, `CA` → `AC`, `N/A` → `--`) and writes
  `sample_data/FINAL.txt` plus an allele count summary `sample_data/STATS.txt`.

---

## Requirements

* [Go](https://go.dev/dl/) 1.15 or newer
* [R](https://cran.r-project.org/) with `Rscript` available on your `PATH`
* R packages: `MASS`, `data.table`, `foreach`, `doParallel`

```r
install.packages(c("MASS", "data.table", "foreach", "doParallel"))
```

Go dependencies (`github.com/gookit/color`, `github.com/lib/pq`) are fetched
automatically by Go modules.

---

## Input files

Put the following files together in one directory (the *project directory*):

| File | Description |
|------|-------------|
| `Project.csv` | GenomeStudio project sample sheet (must contain a `Sample_ID` header line) |
| `submission_template.txt` | Submission template (tab-delimited text export) |
| `FinalReportR.txt` | GenomeStudio Final Report – intensity (R) values (data starts after a line containing `SNP Name`) |
| `FinalReportTheta.txt` | GenomeStudio Final Report – theta values (same format) |

Final Report data lines must be tab-delimited with columns
`SNP Name`, `Sample ID`, `<value>`.

You also need `calls.txt` (included in this repository) in the directory you run
the program from. It maps AB alleles to nucleotide calls per marker.

> These data files are listed in [.gitignore](.gitignore) and are **not** stored
> in the repository.

---

## Running

Clone the repository and run from its root:

```bash
git clone https://github.com/cardinalb/go50k.git
cd go50k

go run main.go -project /path/to/project_directory
```

Or build a binary first:

```bash
go build -o go50k main.go      # go50k.exe on Windows
./go50k -project /path/to/project_directory
```

> Build with `main.go` (or `go build .`); `main.go` is the program's entry point.

| Flag | Description |
|------|-------------|
| `-project` | Directory containing the input files above (no trailing slash) |

The first run on a realistic dataset takes about a minute for the file loading and
up to ~10 minutes for the R clustering, depending on file size.

### Outputs

| File | Location | Description |
|------|----------|-------------|
| `OUT_FinalReportR.txt`, `OUT_FinalReportTheta.txt` | project directory | Transposed marker × sample CSVs passed to R |
| `stats.csv`, `ids.csv`, `parms.csv` | project directory | Per-SNP cluster statistics, per-sample cluster calls and the parameters used |
| `final_calls.txt` | `sample_data/` | ACGT calls |
| `FINAL.txt`, `STATS.txt` | `sample_data/` | Normalised calls and allele frequency counts |

---

## Known limitations

* The post-processing stage reads `sample_data/ids.csv` and writes into
  `sample_data/` using hard-coded relative paths, whereas the R script writes
  `ids.csv` to the project directory. Either use `sample_data` as your project
  directory (`-project ./sample_data`) or copy `ids.csv` into `sample_data/`
  before the conversion step.
* `calls.txt` must be in the current working directory.

---

## Repository layout

```
main.go                  Entry point and AB→ACGT / normalisation functions
cluster50kparallel.R     Parallel k-means genotype clustering
calls.txt                AB → ACGT reference calls per marker
projectfile/             Project.csv loader
submissionfile/          Submission template loader
loaddata/                Final Report loader / transposer
database/                PostgreSQL helpers
setup/                   Terminal UI text and utilities
```

---

## License

Licensed under the Apache License, Version 2.0 – see [LICENSE](LICENSE) and
[NOTICE.txt](NOTICE.txt).
