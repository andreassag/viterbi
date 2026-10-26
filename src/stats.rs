use crate::hmm::ProfileHmm;

pub fn raw_to_bit_score(raw: f32) -> f32 {
    (raw / 0.69314718f32).max(0.0)
}

pub fn compute_evalue(hmm: &ProfileHmm, bit_score: f32, database_size: usize) -> f64 {
    let n = database_size.max(1) as f64;
    let p_val = (-(hmm.lambda * (bit_score - hmm.mu))).exp() as f64;
    (n * p_val).min(1000.0)
}
