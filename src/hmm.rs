use crate::fasta::FastaRecord;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

pub const AMINO_ALPHABET: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";

#[inline(always)]
pub fn aa_to_idx(c: u8) -> usize {
    match c {
        b'A' => 0, b'C' => 1, b'D' => 2, b'E' => 3, b'F' => 4,
        b'G' => 5, b'H' => 6, b'I' => 7, b'K' => 8, b'L' => 9,
        b'M' => 10, b'N' => 11, b'P' => 12, b'Q' => 13, b'R' => 14,
        b'S' => 15, b'T' => 16, b'V' => 17, b'W' => 18, b'Y' => 19,
        _ => 0, // Fallback
    }
}

#[derive(Clone, Debug)]
pub struct ProfileHmm {
    pub name: String,
    pub length: usize, // M nodes
    pub alphabet: Vec<u8>,
    pub match_emissions: Vec<Vec<f32>>,  // [M + 1][alphabet_size]
    pub insert_emissions: Vec<Vec<f32>>, // [M + 1][alphabet_size]
    pub transitions: Vec<[f32; 7]>,      // [M + 1]: MM, MI, MD, IM, II, DM, DD
    pub mu: f32,
    pub lambda: f32,
}

impl ProfileHmm {
    pub fn new(name: &str, length: usize) -> Self {
        let alph = AMINO_ALPHABET.to_vec();
        let k = alph.len();
        let neg_inf = -100.0f32;

        Self {
            name: name.to_string(),
            length,
            alphabet: alph,
            match_emissions: vec![vec![neg_inf; k]; length + 1],
            insert_emissions: vec![vec![neg_inf; k]; length + 1],
            transitions: vec![[neg_inf; 7]; length + 1],
            mu: -10.0,
            lambda: 0.7,
        }
    }

    pub fn build_from_msa(records: &[FastaRecord], name: &str) -> Result<Self, String> {
        if records.is_empty() {
            return Err("Empty MSA".to_string());
        }

        let msa_len = records[0].seq.len();
        for r in records {
            if r.seq.len() != msa_len {
                return Err("MSA sequences have unequal lengths".to_string());
            }
        }

        // Identify match columns (e.g. < 50% gaps)
        let num_seqs = records.len() as f32;
        let mut match_cols = Vec::new();

        for col in 0..msa_len {
            let mut gap_count = 0;
            for r in records {
                if r.seq[col] == b'-' || r.seq[col] == b'.' {
                    gap_count += 1;
                }
            }
            if (gap_count as f32) / num_seqs < 0.5 {
                match_cols.push(col);
            }
        }

        let model_len = match_cols.len();
        let mut hmm = ProfileHmm::new(name, model_len);
        let k = hmm.alphabet.len();
        let pseudo = 0.1f32; // Laplace pseudocount

        // Compute match emissions with Laplace pseudocounts
        for (m_idx, &col) in match_cols.iter().enumerate() {
            let mut counts = vec![pseudo; k];
            let mut total = (k as f32) * pseudo;

            for r in records {
                let b = r.seq[col];
                if b != b'-' && b != b'.' {
                    let idx = aa_to_idx(b);
                    counts[idx] += 1.0;
                    total += 1.0;
                }
            }

            let bg_log = (1.0f32 / (k as f32)).ln();
            for idx in 0..k {
                hmm.match_emissions[m_idx + 1][idx] = (counts[idx] / total).ln() - bg_log;
                hmm.insert_emissions[m_idx + 1][idx] = 0.0;
            }
        }

        // Standard transitions: MM=0.9, MI=0.05, MD=0.05, IM=0.5, II=0.5, DM=0.5, DD=0.5
        for m in 1..=model_len {
            hmm.transitions[m][0] = 0.90f32.ln(); // MM
            hmm.transitions[m][1] = 0.05f32.ln(); // MI
            hmm.transitions[m][2] = 0.05f32.ln(); // MD
            hmm.transitions[m][3] = 0.50f32.ln(); // IM
            hmm.transitions[m][4] = 0.50f32.ln(); // II
            hmm.transitions[m][5] = 0.50f32.ln(); // DM
            hmm.transitions[m][6] = 0.50f32.ln(); // DD
        }

        // Calibration heuristic
        hmm.mu = -((model_len as f32).ln());
        hmm.lambda = 0.693;

        Ok(hmm)
    }

    pub fn save_hmm3(&self, path: &str) -> Result<(), String> {
        let mut file = File::create(path).map_err(|e| format!("Cannot create {}: {}", path, e))?;
        writeln!(file, "HMMER3/f [3.3]").map_err(|e| e.to_string())?;
        writeln!(file, "NAME  {}", self.name).map_err(|e| e.to_string())?;
        writeln!(file, "LENG  {}", self.length).map_err(|e| e.to_string())?;
        writeln!(file, "ALPH  amino").map_err(|e| e.to_string())?;
        writeln!(file, "STATS LOCAL VITERBI   {:.4}  {:.5}", self.mu, self.lambda).map_err(|e| e.to_string())?;
        writeln!(file, "HMM").map_err(|e| e.to_string())?;
        writeln!(file, "        A        C        D        E        F        G        H        I        K        L        M        N        P        Q        R        S        T        V        W        Y").map_err(|e| e.to_string())?;
        writeln!(file, "      m->m     m->i     m->d     i->m     i->i     d->m     d->d").map_err(|e| e.to_string())?;

        for m in 1..=self.length {
            write!(file, "  {:5d}", m).map_err(|e| e.to_string())?;
            for val in &self.match_emissions[m] {
                write!(file, " {:7.2}", val).map_err(|e| e.to_string())?;
            }
            writeln!(file).map_err(|e| e.to_string())?;

            write!(file, "       ").map_err(|e| e.to_string())?;
            for val in &self.insert_emissions[m] {
                write!(file, " {:7.2}", val).map_err(|e| e.to_string())?;
            }
            writeln!(file).map_err(|e| e.to_string())?;

            write!(file, "       ").map_err(|e| e.to_string())?;
            for val in &self.transitions[m] {
                write!(file, " {:7.2}", val).map_err(|e| e.to_string())?;
            }
            writeln!(file).map_err(|e| e.to_string())?;
        }

        writeln!(file, "//").map_err(|e| e.to_string())?;
        Ok(())
    }

    pub fn load_hmm3(path: &str) -> Result<Self, String> {
        let file = File::open(path).map_err(|e| format!("Cannot open {}: {}", path, e))?;
        let reader = BufReader::new(file);

        let mut name = "unnamed".to_string();
        let mut length = 0;
        let mut mu = -10.0f32;
        let mut lambda = 0.7f32;

        let mut lines = reader.lines();
        while let Some(Ok(line)) = lines.next() {
            let line = line.trim();
            if line.starts_with("NAME") {
                name = line[4..].trim().to_string();
            } else if line.starts_with("LENG") {
                length = line[4..].trim().parse().unwrap_or(0);
            } else if line.starts_with("STATS LOCAL VITERBI") {
                let parts: Vec<&str> = line[19..].split_whitespace().collect();
                if parts.len() >= 2 {
                    mu = parts[0].parse().unwrap_or(-10.0);
                    lambda = parts[1].parse().unwrap_or(0.7);
                }
            } else if line == "HMM" {
                break;
            }
        }

        if length == 0 {
            return Err("Invalid HMM file: missing or zero LENG".to_string());
        }

        let mut hmm = ProfileHmm::new(&name, length);
        hmm.mu = mu;
        hmm.lambda = lambda;

        // Skip 2 header lines (alphabet and transition headers)
        lines.next();
        lines.next();

        let mut cur_node = 1;
        while cur_node <= length {
            // Match emission line
            if let Some(Ok(line)) = lines.next() {
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() >= 21 {
                    for i in 0..20 {
                        hmm.match_emissions[cur_node][i] = parts[i + 1].parse().unwrap_or(-100.0);
                    }
                }
            }
            // Insert emission line
            if let Some(Ok(line)) = lines.next() {
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() >= 20 {
                    for i in 0..20 {
                        hmm.insert_emissions[cur_node][i] = parts[i].parse().unwrap_or(-100.0);
                    }
                }
            }
            // Transition line
            if let Some(Ok(line)) = lines.next() {
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() >= 7 {
                    for i in 0..7 {
                        hmm.transitions[cur_node][i] = parts[i].parse().unwrap_or(-100.0);
                    }
                }
            }
            cur_node += 1;
        }

        Ok(hmm)
    }
}
