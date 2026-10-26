#!/usr/bin/env bash
set -euo pipefail

echo "=== 1. Verifying Profile HMM Discrimination Logic ==="
python3 scripts/verify_logic.py

if command -v cargo >/dev/null 2>&1; then
    echo -e "\n=== 2. Running Rust Viterbi CLI ==="
    cargo run --release -- build --msa data/sample_msa.fa --output kinase.hmm
    cargo run --release -- search --hmm kinase.hmm --targets data/sample_targets.fa --domtblout hits.tsv --show-alignment
    rm -f kinase.hmm hits.tsv
else
    echo -e "\n[Note: cargo toolchain not yet installed; algorithmic verification ran via Python suite]"
fi

echo -e "\nExample run completed successfully!"
