#!/usr/bin/env bash
set -euo pipefail

echo "=== Benchmark: Profile HMM Viterbi Search Throughput ==="

python3 -c '
import time, random
from verify_logic import build_simple_hmm, viterbi_score, AMINO_ALPHABET

msa = []
for _ in range(10):
    msa.append("".join(random.choice(AMINO_ALPHABET) for _ in range(100)))

hmm = build_simple_hmm(msa)

for target_len in [200, 500, 1000]:
    target = "".join(random.choice(AMINO_ALPHABET) for _ in range(target_len))

    t0 = time.perf_counter()
    for _ in range(20):
        viterbi_score(hmm, target)
    elapsed = (time.perf_counter() - t0) / 20

    print(f"Model 100 nodes x Target {target_len:4d} aa: Viterbi Latency = {elapsed*1000:6.2f} ms")
'

echo "Benchmark completed."
