use crate::hmm::{aa_to_idx, ProfileHmm};

#[derive(Clone, Debug)]
pub struct ViterbiResult {
    pub raw_score: f32,
    pub bit_score: f32,
    pub evalue: f64,
    pub query_start: usize,
    pub query_end: usize,
    pub model_start: usize,
    pub model_end: usize,
    pub state_path: Vec<char>, // 'M', 'I', 'D'
}

pub fn viterbi_scalar(hmm: &ProfileHmm, query: &[u8]) -> ViterbiResult {
    let l = query.len();
    let m = hmm.length;
    let neg_inf = -1000.0f32;

    if l == 0 || m == 0 {
        return ViterbiResult {
            raw_score: neg_inf,
            bit_score: 0.0,
            evalue: 1000.0,
            query_start: 0,
            query_end: 0,
            model_start: 0,
            model_end: 0,
            state_path: Vec::new(),
        };
    }

    // DP matrices [M + 1][L + 1]
    let mut dp_m = vec![vec![neg_inf; l + 1]; m + 1];
    let mut dp_i = vec![vec![neg_inf; l + 1]; m + 1];
    let mut dp_d = vec![vec![neg_inf; l + 1]; m + 1];

    dp_m[0][0] = 0.0;

    let mut best_score = neg_inf;
    let mut best_k = 0;
    let mut best_j = 0;

    for j in 1..=l {
        let sym = aa_to_idx(query[j - 1]);

        for k in 1..=m {
            // Match state
            let t = hmm.transitions[k];
            let em = hmm.match_emissions[k][sym];

            // Local start possibility at any model position
            let from_begin = if k == 1 { 0.0 } else { neg_inf };

            let m_prev = dp_m[k - 1][j - 1] + t[0];
            let i_prev = dp_i[k - 1][j - 1] + t[3];
            let d_prev = dp_d[k - 1][j - 1] + t[5];
            let best_in = m_prev.max(i_prev).max(d_prev).max(from_begin);

            dp_m[k][j] = best_in + em;

            // Insert state
            let ei = hmm.insert_emissions[k][sym];
            let m_to_i = dp_m[k][j - 1] + t[1];
            let i_to_i = dp_i[k][j - 1] + t[4];
            dp_i[k][j] = m_to_i.max(i_to_i) + ei;

            // Delete state
            let m_to_d = dp_m[k - 1][j] + t[2];
            let d_to_d = dp_d[k - 1][j] + t[6];
            dp_d[k][j] = m_to_d.max(d_to_d);

            let cur_max = dp_m[k][j].max(dp_i[k][j]).max(dp_d[k][j]);
            if cur_max > best_score {
                best_score = cur_max;
                best_k = k;
                best_j = j;
            }
        }
    }

    // Backtrace
    let mut cur_k = best_k;
    let mut cur_j = best_j;
    let mut cur_state = if dp_m[best_k][best_j] >= dp_i[best_k][best_j] && dp_m[best_k][best_j] >= dp_d[best_k][best_j] {
        'M'
    } else if dp_i[best_k][best_j] >= dp_d[best_k][best_j] {
        'I'
    } else {
        'D'
    };

    let mut path = Vec::new();

    while cur_k > 0 && cur_j > 0 {
        if cur_state == 'M' {
            path.push('M');
            let t = hmm.transitions[cur_k];
            let from_begin = if cur_k == 1 { 0.0 } else { neg_inf };

            let m_prev = dp_m[cur_k - 1][cur_j - 1] + t[0];
            let i_prev = dp_i[cur_k - 1][cur_j - 1] + t[3];
            let d_prev = dp_d[cur_k - 1][cur_j - 1] + t[5];

            if from_begin >= m_prev && from_begin >= i_prev && from_begin >= d_prev {
                cur_k -= 1;
                cur_j -= 1;
                break;
            } else if m_prev >= i_prev && m_prev >= d_prev {
                cur_state = 'M';
            } else if i_prev >= d_prev {
                cur_state = 'I';
            } else {
                cur_state = 'D';
            }
            cur_k -= 1;
            cur_j -= 1;
        } else if cur_state == 'I' {
            path.push('I');
            let t = hmm.transitions[cur_k];
            let m_to_i = dp_m[cur_k][cur_j - 1] + t[1];
            let i_to_i = dp_i[cur_k][cur_j - 1] + t[4];
            if m_to_i >= i_to_i {
                cur_state = 'M';
            } else {
                cur_state = 'I';
            }
            cur_j -= 1;
        } else {
            path.push('D');
            let t = hmm.transitions[cur_k];
            let m_to_d = dp_m[cur_k - 1][cur_j] + t[2];
            let d_to_d = dp_d[cur_k - 1][cur_j] + t[6];
            if m_to_d >= d_to_d {
                cur_state = 'M';
            } else {
                cur_state = 'D';
            }
            cur_k -= 1;
        }
    }

    path.reverse();
    let q_start = cur_j;
    let m_start = cur_k;

    // Bit-score conversion
    let bit_score = (best_score / 0.693147).max(0.0);
    // E-value calculation via Gumbel distribution
    let evalue = (-(hmm.lambda * (bit_score - hmm.mu))).exp() as f64;

    ViterbiResult {
        raw_score: best_score,
        bit_score,
        evalue: evalue.min(1000.0),
        query_start: q_start,
        query_end: best_j,
        model_start: m_start,
        model_end: best_k,
        state_path: path,
    }
}
