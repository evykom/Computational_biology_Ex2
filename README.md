# Computational Biology Exercise 2

This project optimizes a DNA sequence for human insulin so it can be expressed efficiently in *E. coli*. The amino-acid sequence stays fixed; the genetic algorithm chooses synonymous codons and balances codon usage, GC content, and repeated codons.

## Requirements

- Python 3.10 or newer
- No external packages

## How to run

Clone the repository and move into the project root before running anything. All commands below must be run from the root directory, where `main.py`, `genetic_code.json`, `ecoli_frequencies.json`, and `insulin_target.txt` are located.

```bash
git https://github.com/evykom/Computational_biology_Ex2.git
cd Computational_biology_Ex2
```

Basic genetic algorithm:

```bash
python main.py
```

Run a specific strategy:

```bash
python main.py --strategy basic
python main.py --strategy baldwin
python main.py --strategy lamarckian
```

Run all three strategies in one command:

```bash
python main.py --strategy all
```

Save per-generation results to CSV files:

```bash
python main.py --strategy all --csv
```

Change the local-search budget for Baldwin and Lamarckian runs:

```bash
python main.py --strategy baldwin --local-search-budget 3
```

## Command-line options


| Flag                    | Description                                | Default |
| ----------------------- | ------------------------------------------ | ------- |
| `--strategy`            | `basic`, `baldwin`, `lamarckian`, or `all` | `basic` |
| `--local-search-budget` | Number of lifetime local-search fixes      | `5`     |
| `--csv`                 | Write per-generation results to CSV        | off     |


## Output

Each run prints the best DNA sequence, its fitness score, the total number of fitness evaluations, and the best score for every generation.

With `--csv`, the program creates a `csv_files` folder in the project root and saves one file per strategy there, for example `csv_files/darwinian_basic_generations.csv`. Each file contains the generation number, average population fitness, best fitness in that generation, and the best fitness found so far.