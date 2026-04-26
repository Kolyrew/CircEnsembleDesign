# CircEnsembleDesign: Design of Circular mRNAs Based on Minimization of the Free Energy of Ensemble Secondary Structures

CircEnsembleDesign is designed for circular mRNA sequence design from protein input, with optimization of ensemble free energy and additional control of undesirable base pairing between IRES and the rest of the molecule.
The algorithm keeps the core probabilistic lattice parsing and gradient-based optimization workflow from EnsembleDesign, and extends it with circular-RNA-specific constraints: fixed IRES prefix handling, IRES-vs-rest cross-pair penalty in DP scoring, and final circular structure filtering using circular RNAfold.


## Dependencies
- Clang 11.0.0 (or above) or GCC 4.8.5 (or above)
- Python 3
- ViennaRNA Python bindings (`import RNA`) for `RNAfold` wrapper used in circular post-filtering

## To Compile
Run the following commands:

```
chmod +x setup.sh
./setup.sh
```

## Usage

The mRNA Design program is designed to process protein sequences provided in a FASTA file and apply our mRNA Design algorithm to each sequence. The results will be stored in `output_dir`, with each sequence having its own subfolder containing the raw outputs from all executions. The optimal solution for each sequence, determined across all runs, will be displayed to the user.

### Command Syntax

Use the following command format to run the mRNA Design program:

```
python EnsembleDesign.py [--fasta <path>] [--output_dir <path>] [--beam_size <int>] [--lr <float>] [--epsilon <float>] [--num_iters <int>] [--num_runs <int>] [--num_threads <int>] [--ires <string>] [--ires_orf_lambda <float>] [--max_cross_pairs <int>] [--cross_pair_prob_threshold <float>]
```

- `--fasta <path>`: Specifies the path to the input protein FASTA file. The default is `./examples.fasta`.
- `--output_dir <path>`: Sets the directory for saving output files. The default is `./outputs`.
- `--beam_size <int>`: Determines the beam size for beam pruning. The default is `200`. A smaller beam size can speed up the process but may result in search errors.
- `--lr <float>`: Sets the learning rate for projected gradient descent. The default is `0.03`.
- `--epsilon <float>`: Defines the epsilon parameter for soft-MFE initialization. The default is `0.5`.
- `--num_iters <int>`: Specifies the number of optimization steps. The default is `30`.
- `--num_runs <int>`: Indicates the number of execution runs per sample file. The default is `20`. More runs increase the likelihood of finding optimal solutions but require more computational resources.
- `--num_threads <int>`: Sets the number of threads in the thread pool. The default is `8`. Adjust this based on your CPU's core count to optimize parallel processing without overloading your system.
- `--ires <string>`: Fixed RNA prefix (IRES) prepended to the designed sequence.
- `--ires_orf_lambda <float>`: Soft penalty factor for IRES-vs-rest pairing inside DP (`0 < lambda <= 1`; default `1.0` means no penalty). Lower values enforce stronger discouragement.
- `--max_cross_pairs <int>`: Maximum allowed number of IRES-vs-rest cross-pairs in final circular structure filtering.
- `--cross_pair_prob_threshold <float>`: Kept for compatibility with previous probabilistic checks. Current final filtering uses circular `RNAfold` structure.

### Current Design + Filtering Pipeline

For each protein sequence:

1. `num_runs` independent optimization runs are executed in C++ (`bin/EnsembleDesign`).
2. Each run applies the ensemble objective and (optionally) IRES penalty in DP via `--ires_orf_lambda`.
3. Final candidate sequences are re-evaluated in Python and folded with circular RNA mode (`python3 ./RNAfold --circ --noPS`).
4. Cross-pairs are counted as base pairs with one index in IRES and the other outside IRES.
5. Selection rule:
   - Prefer the lowest EFE candidate among those with `cross_pairs <= max_cross_pairs`.
   - If none satisfy the threshold, choose the best available candidate with minimum cross-pairs (and best EFE tie-break).

### Example Command and Expected Output

To run the program with a specific set of parameters, you can use a command similar to the following (note that we are using smaller parameters to make it run faster):

```
python EnsembleDesign.py \
  --fasta ./data/examples.fasta \
  --output_dir ./outputs \
  --beam_size 100 \
  --lr 0.03 \
  --epsilon 0.5 \
  --num_iters 10 \
  --num_runs 2 \
  --num_threads 2 \
  --ires "UUCGUUGCUUUUUGUAGUAUAAUUAAAUAUUUGUCAUAUAAGAGAUUGGUUAGAGAUUUGUUCUUUGUUUGA" \
  --ires_orf_lambda 0.1 \
  --max_cross_pairs 5
```

The expected output format is:

```
>seq1|Ensemble Free Energy: -11.52 kcal/mol|Cross-pairs: 5
UUCGUUGCUUUUUGUAGUAUAAAUGAAUACUUACCAUAUUACUUUGCCUUGGCCUCCUAGCAAUAAUAGGUAC
```

Per-run logs are written into `output_dir/<sequence_id>/run_*.txt`, and final filtering details (including circular structure and `OK/REJECTED`) are written to `output_dir/<sequence_id>/results.txt`.

Note that following LinearDesign, we use `-V -d0 -b0` options in LinearPartition to evaluate ensemble free energy, which uses the Vienna energy model with no dangling ends and no beam search (although the mRNA design code itself uses beam search).

## Other Files

Alongside the main application, this repository includes additional files that are useful for testing and understanding the capabilities of the mRNA Design tool:

- `data/examples.fasta`: This file contains 3 example protein sequences. It is designed to offer a quick and straightforward way to test the functionality of the software with pre-defined input. This file serves as the default input for the tool.

- `data/uniprot.fasta`: Contains 20 protein sequences from the UniProt database used in our experiments.

- `data/covid_spike.fasta`: Contains the SARS-CoV-2 spike protein sequence used in our experiments, taken from the LinearDesign paper (Zhang et al., Nature 2023).

- `coding_wheel.txt`: RNA codon table (i.e., the genetic code).

- `codon_usage_freq_table_human.csv`: A human codon usage frequency table from the LinearDesign paper, which is used by LinearDesign code to compute the codon adapation index (CAI). It is originally from [GenScript](https://www.genscript.com/tools/codon-frequency-table) (choose human instead of yeast).

- `tools/LinearDesign`: [LinearDesign codebase](https://github.com/LinearDesignSoftware/LinearDesign) (baseline, for min MFE design).

- `tools/LinearPartition`: [LinearPartition codebase](https://github.com/LinearFold/LinearPartition) (to calculate ensemble free energy). Note that we use `-V -d0 -b0` options to evaluate energies, following LinearDesign.