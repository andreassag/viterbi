use std::fs::File;
use std::io::{BufRead, BufReader};

#[derive(Clone, Debug)]
pub struct FastaRecord {
    pub name: String,
    pub seq: Vec<u8>,
}

pub fn read_fasta(path: &str) -> Result<Vec<FastaRecord>, String> {
    let file = File::open(path).map_err(|e| format!("Cannot open {}: {}", path, e))?;
    let reader = BufReader::new(file);

    let mut records = Vec::new();
    let mut cur_name = String::new();
    let mut cur_seq = Vec::new();

    for line_res in reader.lines() {
        let line = line_res.map_err(|e| e.to_string())?;
        let line = line.trim();
        if line.is_empty() {
            continue;
        }

        if line.starts_with('>') {
            if !cur_name.is_empty() {
                records.push(FastaRecord {
                    name: cur_name,
                    seq: cur_seq,
                });
                cur_seq = Vec::new();
            }
            cur_name = line[1..].split_whitespace().next().unwrap_or("seq").to_string();
        } else {
            for b in line.bytes() {
                if b > b' ' {
                    cur_seq.push(b.to_ascii_uppercase());
                }
            }
        }
    }

    if !cur_name.is_empty() {
        records.push(FastaRecord {
            name: cur_name,
            seq: cur_seq,
        });
    }

    Ok(records)
}
