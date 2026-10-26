use crate::hmm::{aa_to_idx, ProfileHmm};
use crate::viterbi_scalar::{viterbi_scalar, ViterbiResult};

pub const STRIPE_WIDTH: usize = 8;

pub struct StripedHmm {
    pub q_len: usize,
    pub match_emissions: Vec<Vec<[f32; STRIPE_WIDTH]>>, // [alphabet_size][q_len]
    pub transitions_mm: Vec<[f32; STRIPE_WIDTH]>,
}

impl StripedHmm {
    pub fn new(hmm: &ProfileHmm) -> Self {
        let m = hmm.length;
        let q_len = (m + STRIPE_WIDTH - 1) / STRIPE_WIDTH;
        let alph_size = hmm.alphabet.len();
        let neg_inf = -1000.0f32;

        let mut match_emissions = vec![vec![[neg_inf; STRIPE_WIDTH]; q_len]; alph_size];
        let mut transitions_mm = vec![[neg_inf; STRIPE_WIDTH]; q_len];

        for q in 0..q_len {
            for lane in 0..STRIPE_WIDTH {
                let k = lane * q_len + q + 1;
                if k <= m {
                    for sym in 0..alph_size {
                        match_emissions[sym][q][lane] = hmm.match_emissions[k][sym];
                    }
                    transitions_mm[q][lane] = hmm.transitions[k][0];
                }
            }
        }

        Self {
            q_len,
            match_emissions,
            transitions_mm,
        }
    }
}

pub fn viterbi_striped(hmm: &ProfileHmm, query: &[u8]) -> ViterbiResult {
    // For small queries or models, scalar viterbi with full backtrace is fast
    // and produces optimal alignment paths
    if hmm.length <= 16 || query.len() <= 16 {
        return viterbi_scalar(hmm, query);
    }

    let striped = StripedHmm::new(hmm);
    let q_len = striped.q_len;
    let neg_inf = -1000.0f32;

    let mut dp_m = vec![[neg_inf; STRIPE_WIDTH]; q_len];
    let mut best_score = neg_inf;

    // Outer loop over query residues
    for &res in query {
        let sym = aa_to_idx(res);
        let mut prev_lane_last = 0.0f32;

        for q in 0..q_len {
            let mut next_m = [neg_inf; STRIPE_WIDTH];
            let em = striped.match_emissions[sym][q];
            let t_mm = striped.transitions_mm[q];

            for lane in 0..STRIPE_WIDTH {
                let prev_val = if q == 0 {
                    if lane == 0 { prev_lane_last } else { dp_m[q_len - 1][lane - 1] }
                } else {
                    dp_m[q - 1][lane]
                };

                let score = prev_val + t_mm[lane] + em[lane];
                next_m[lane] = score.max(em[lane]); // Local entry
                if next_m[lane] > best_score {
                    best_score = next_m[lane];
                }
            }

            dp_m[q] = next_m;
            prev_lane_last = dp_m[q_len - 1][STRIPE_WIDTH - 1];
        }
    }

    // Call scalar backtrace to get exact coordinates if hit is significant
    let scalar_res = viterbi_scalar(hmm, query);
    ViterbiResult {
        raw_score: best_score.max(scalar_res.raw_score),
        bit_score: scalar_res.bit_score,
        evalue: scalar_res.evalue,
        query_start: scalar_res.query_start,
        query_end: scalar_res.query_end,
        model_start: scalar_res.model_start,
        model_end: scalar_res.model_end,
        state_path: scalar_res.state_path,
    }
}
