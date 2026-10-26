pub mod fasta;
pub mod math;
pub mod hmm;
pub mod viterbi_scalar;
pub mod viterbi_striped;
pub mod forward;
pub mod stats;

pub use fasta::{read_fasta, FastaRecord};
pub use math::log_sum_exp;
pub use hmm::{ProfileHmm, AMINO_ALPHABET, aa_to_idx};
pub use viterbi_scalar::{viterbi_scalar, ViterbiResult};
pub use viterbi_striped::viterbi_striped;
pub use forward::forward_score;
pub use stats::{raw_to_bit_score, compute_evalue};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_msa_build_and_viterbi() {
        let msa = vec![
            FastaRecord { name: "s1".to_string(), seq: b"ACDEFGHIKLMN".to_vec() },
            FastaRecord { name: "s2".to_string(), seq: b"ACDEFGHIKLMN".to_vec() },
            FastaRecord { name: "s3".to_string(), seq: b"ACDEFGHIKLMN".to_vec() },
        ];
        let hmm = ProfileHmm::build_from_msa(&msa, "test_model").unwrap();
        assert_eq!(hmm.length, 12);

        // Identical query should score high
        let res = viterbi_scalar(&hmm, b"ACDEFGHIKLMN");
        assert!(res.raw_score > -20.0);
        assert!(res.bit_score > 0.0);
        assert_eq!(res.query_start, 0);
        assert_eq!(res.query_end, 12);
    }

    #[test]
    fn test_striped_vs_scalar() {
        let msa = vec![
            FastaRecord { name: "s1".to_string(), seq: b"ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY".to_vec() },
            FastaRecord { name: "s2".to_string(), seq: b"ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY".to_vec() },
        ];
        let hmm = ProfileHmm::build_from_msa(&msa, "long_model").unwrap();

        let target = b"ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY";
        let res_scalar = viterbi_scalar(&hmm, target);
        let res_striped = viterbi_striped(&hmm, target);

        assert!((res_scalar.raw_score - res_striped.raw_score).abs() < 5.0);
    }
}
