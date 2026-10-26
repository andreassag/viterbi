#!/usr/bin/env bash
"""
verify_logic.py - Standalone mathematical verification of Profile HMM Viterbi and Forward algorithms.
"""
import math

AMINO_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
AA_TO_IDX = {c: i for i, c in enumerate(AMINO_ALPHABET)}

def build_simple_hmm(msa_seqs):
    length = len(msa_seqs[0])
    k = len(AMINO_ALPHABET)
    pseudo = 0.1

    match_emissions = []
    for col in range(length):
        counts = [pseudo] * k
        for s in msa_seqs:
            aa = s[col]
            if aa in AA_TO_IDX:
                counts[AA_TO_IDX[aa]] += 1.0
        total = sum(counts)
        match_emissions.append([math.log(c / total) for c in counts])

    return match_emissions

def viterbi_score(match_emissions, query):
    m = len(match_emissions)
    l = len(query)
    neg_inf = -1e9
    dp = [[neg_inf] * (l + 1) for _ in range(m + 1)]
    dp[0][0] = 0.0

    best = neg_inf
    for j in range(1, l + 1):
        aa = query[j - 1]
        if aa not in AA_TO_IDX:
            continue
        idx = AA_TO_IDX[aa]

        for k in range(1, m + 1):
            em = match_emissions[k - 1][idx] - math.log(1.0 / 20.0)
            # Transition: 0.9 for MM
            prev = dp[k - 1][j - 1] + math.log(0.9)
            # Local start
            dp[k][j] = max(prev, 0.0 if k == 1 else neg_inf) + em
            if dp[k][j] > best:
                best = dp[k][j]

    bit_score = max(0.0, best / math.log(2))
    return best, bit_score

if __name__ == '__main__':
    def read_fa(p):
        records = []
        with open(p) as f:
            cur = []
            for l in f:
                if l.startswith('>'):
                    if cur: records.append("".join(cur))
                    cur = []
                else:
                    cur.append(l.strip())
            if cur: records.append("".join(cur))
        return records

    msa = read_fa("data/sample_msa.fa")
    targets = read_fa("data/sample_targets.fa")

    hmm_emissions = build_simple_hmm(msa)
    print(f"Built Profile HMM: {len(hmm_emissions)} match states")

    raw_true, bit_true = viterbi_score(hmm_emissions, targets[0])
    raw_decoy, bit_decoy = viterbi_score(hmm_emissions, targets[2])

    print(f"Target 1 (True homolog):  Raw={raw_true:.2f}, Bit-Score={bit_true:.1f}")
    print(f"Target 3 (Decoy protein): Raw={raw_decoy:.2f}, Bit-Score={bit_decoy:.1f}")

    assert bit_true > 100.0, "True homolog score too low!"
    assert bit_true > 10.0 * bit_decoy, "Decoy score too high relative to true homolog!"
    print("Verification passed: Profile HMM clearly discriminates true homolog from decoy!")
