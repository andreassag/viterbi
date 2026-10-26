use crate::hmm::{aa_to_idx, ProfileHmm};
use crate::math::log_sum_exp;

pub fn forward_score(hmm: &ProfileHmm, query: &[u8]) -> f32 {
    let l = query.len();
    let m = hmm.length;
    let neg_inf = -1000.0f32;

    if l == 0 || m == 0 {
        return neg_inf;
    }

    let mut dp_m = vec![vec![neg_inf; l + 1]; m + 1];
    let mut dp_i = vec![vec![neg_inf; l + 1]; m + 1];
    let mut dp_d = vec![vec![neg_inf; l + 1]; m + 1];

    dp_m[0][0] = 0.0;

    for j in 1..=l {
        let sym = aa_to_idx(query[j - 1]);

        for k in 1..=m {
            let t = hmm.transitions[k];
            let em = hmm.match_emissions[k][sym];

            // Match state sum
            let m_prev = dp_m[k - 1][j - 1] + t[0];
            let i_prev = dp_i[k - 1][j - 1] + t[3];
            let d_prev = dp_d[k - 1][j - 1] + t[5];
            let sum_in = log_sum_exp(log_sum_exp(m_prev, i_prev), d_prev);
            dp_m[k][j] = sum_in + em;

            // Insert state sum
            let ei = hmm.insert_emissions[k][sym];
            let m_to_i = dp_m[k][j - 1] + t[1];
            let i_to_i = dp_i[k][j - 1] + t[4];
            dp_i[k][j] = log_sum_exp(m_to_i, i_to_i) + ei;

            // Delete state sum
            let m_to_d = dp_m[k - 1][j] + t[2];
            let d_to_d = dp_d[k - 1][j] + t[6];
            dp_d[k][j] = log_sum_exp(m_to_d, d_to_d);
        }
    }

    // Sum over all terminal match states
    let mut total = neg_inf;
    for k in 1..=m {
        total = log_sum_exp(total, dp_m[k][l]);
    }
    total
}
