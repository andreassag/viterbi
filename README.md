> [!IMPORTANT]
> **Archival Notice**
> This repository is archived and is no longer actively maintained.

# viterbi: Striped Profile HMM Search Engine

[![Rust](https://img.shields.io/badge/rust-2021-orange.svg)](https://www.rust-lang.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Zero Dependencies](https://img.shields.io/badge/dependencies-none-brightgreen.svg)]()

`viterbi` is a high-performance profile Hidden Markov Model (HMM) search and alignment engine written in pure Rust. It implements HMMER3-compatible profile modeling, fast log-space arithmetic with a 4,096-entry $\ln(1 + e^{-\Delta})$ lookup table, scalar local/glocal Viterbi search with backtrace, Farrar-style striped vector processing, and extreme-value Gumbel calibration for bit-scores and E-values.

---

## Key Features

1. **Precomputed Branchless $\ln(1 + e^{-\Delta})$ Table**:
   - 4,096-entry lookup table for high-speed `log_sum_exp(a, b)` evaluations, enabling fast Forward scoring without calling `f32::exp` or `f32::ln` in the inner loop.
2. **Dual DP Engines (Scalar & Striped)**:
   - **Scalar Viterbi**: Full local profile HMM recurrence with state backtrace (`M`, `I`, `D`), yielding domain boundaries and state alignments.
   - **Farrar-Style Striped Viterbi**: Parallelized model states into striped lanes of width $W = 8$, boosting cache locality and SIMD auto-vectorization across query sequences.
3. **HMMER3 ASCII Format (`.hmm`) & MSA Builder**:
   - Parses standard ASCII HMMER3/f model files (`NAME`, `LENG`, `ALPH`, `STATS`, match/insert/transition matrices).
   - Builds calibrated profile HMMs directly from multi-sequence alignment (MSA) FASTA files with Laplace pseudocounts.
4. **Statistical Significance**:
   - Computes log-odds **bit-scores** and statistical **E-values** calibrated via Gumbel extreme value distribution parameters ($\mu, \lambda$).
   - Generates tabular domain hit tables (`--domtblout`).

---

## Architectural Overview

```
Multiple Sequence Alignment (MSA)
              |
       [ viterbi build ]
              |
     Laplace Pseudocounts + Background Normalization
              |
       Profile HMM (.hmm)
              |
              v
       [ Target Sequences FASTA ]
              |
       +------+------+
       |             |
 [ Striped DP ] [ Scalar DP ] (W=8 vectorization & backtrace)
       |             |
       +------+------+
              |
   [ Gumbel Extreme Value Statistics ]
   - Bit-Score S = (Raw / ln 2)
   - E-value E = N * exp(-lambda * (S - mu))
              |
      +-------+-------+
      |               |
Standard Output   Tabular Hits (domtblout)
```

---

## Build & Test

### Building via Cargo
```bash
cargo build --release
cargo test --release
```

### Running via Makefile
```bash
make test       # Runs unit tests or standalone verification suite
make all        # Builds release binary
```

---

## CLI Usage

### 1. Build Profile HMM from MSA FASTA
```bash
cargo run --release -- build --msa data/sample_msa.fa --output kinase.hmm
```

### 2. Search Profile HMM against Protein Targets
```bash
cargo run --release -- search --hmm kinase.hmm --targets data/sample_targets.fa --domtblout hits.tsv --show-alignment
```

Output:
```text
Target                      Length      Score      E-value hmm_from   hmm_to ali_from   ali_to
------------------------------------------------------------------------------------------
target_kinase_true              59      152.1     1.42e-44        0       41        4       45
  States: MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
target_kinase_variant           50       84.3     2.15e-23        0       41        4       45
  States: MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
```

---

## Verification & Benchmarks

Run automated verification:
```bash
./scripts/test_example.sh
./scripts/benchmark_compare.sh
```

### Performance on Model (100 nodes):
| Target Length | Viterbi Latency |
| :--- | :--- |
| **200 aa** | **4.88 ms** |
| **500 aa** | **12.76 ms** |
| **1,000 aa** | **26.74 ms** |

---

## License

MIT License. Copyright (c) 2026 Andreas Sagen.
