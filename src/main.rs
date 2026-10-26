use std::env;
use std::fs::File;
use std::io::Write;
use std::process;
use viterbi::{
    compute_evalue, read_fasta, viterbi_scalar, viterbi_striped, ProfileHmm,
};

fn print_usage() {
    eprintln!(
        "Usage: viterbi <command> [options]\n\n\
         Commands:\n\
           build --msa <msa.fa> --output <model.hmm>\n\
                 Build Profile HMM from multiple sequence alignment FASTA\n\
           search --hmm <model.hmm> --targets <targets.fa> [options]\n\
                 Search Profile HMM against target sequences\n\n\
         Search Options:\n\
           --mode <striped|scalar>   Search algorithm mode (default: striped)\n\
           --evalue <float>          E-value reporting threshold (default: 10.0)\n\
           --domtblout <file>        Output tabular domain hits TSV file\n\
           --show-alignment          Print alignment coordinates\n"
    );
}

fn cmd_build(args: &[String]) {
    let mut msa_path = String::new();
    let mut out_path = String::new();

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--msa" => {
                i += 1;
                if i < args.len() { msa_path = args[i].clone(); }
            }
            "--output" | "-o" => {
                i += 1;
                if i < args.len() { out_path = args[i].clone(); }
            }
            _ => {}
        }
        i += 1;
    }

    if msa_path.is_empty() || out_path.is_empty() {
        eprintln!("Error: --msa and --output are required for build.");
        process::exit(1);
    }

    let records = read_fasta(&msa_path).unwrap_or_else(|e| {
        eprintln!("Error reading MSA: {}", e);
        process::exit(1);
    });

    println!("[viterbi] Building HMM from {} sequences in {}...", records.len(), msa_path);
    let hmm = ProfileHmm::build_from_msa(&records, "msa_model").unwrap_or_else(|e| {
        eprintln!("Error building HMM: {}", e);
        process::exit(1);
    });

    hmm.save_hmm3(&out_path).unwrap_or_else(|e| {
        eprintln!("Error saving HMM: {}", e);
        process::exit(1);
    });

    println!("[viterbi] Saved Profile HMM (length {} nodes) to {}", hmm.length, out_path);
}

fn cmd_search(args: &[String]) {
    let mut hmm_path = String::new();
    let mut targets_path = String::new();
    let mut mode = "striped".to_string();
    let mut max_evalue = 10.0f64;
    let mut domtblout_path = String::new();
    let mut show_alignment = false;

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--hmm" => {
                i += 1;
                if i < args.len() { hmm_path = args[i].clone(); }
            }
            "--targets" | "-t" => {
                i += 1;
                if i < args.len() { targets_path = args[i].clone(); }
            }
            "--mode" => {
                i += 1;
                if i < args.len() { mode = args[i].clone(); }
            }
            "--evalue" => {
                i += 1;
                if i < args.len() { max_evalue = args[i].parse().unwrap_or(10.0); }
            }
            "--domtblout" => {
                i += 1;
                if i < args.len() { domtblout_path = args[i].clone(); }
            }
            "--show-alignment" => {
                show_alignment = true;
            }
            _ => {}
        }
        i += 1;
    }

    if hmm_path.is_empty() || targets_path.is_empty() {
        eprintln!("Error: --hmm and --targets are required for search.");
        process::exit(1);
    }

    let hmm = ProfileHmm::load_hmm3(&hmm_path).unwrap_or_else(|e| {
        eprintln!("Error loading HMM: {}", e);
        process::exit(1);
    });

    let targets = read_fasta(&targets_path).unwrap_or_else(|e| {
        eprintln!("Error reading targets: {}", e);
        process::exit(1);
    });

    println!("[viterbi] Loaded HMM '{}' ({} nodes), searching {} targets...", hmm.name, hmm.length, targets.len());

    let mut domtbl_file = if !domtblout_path.is_empty() {
        let mut f = File::create(&domtblout_path).unwrap();
        writeln!(f, "# target_name\ttarget_len\tquery_name\tquery_len\tbit_score\tevalue\thmm_from\thmm_to\taln_from\taln_to").unwrap();
        Some(f)
    } else {
        None
    };

    println!("\n{:20} {:>10} {:>10} {:>12} {:>8} {:>8} {:>8} {:>8}", "Target", "Length", "Score", "E-value", "hmm_from", "hmm_to", "ali_from", "ali_to");
    println!("{:-<90}", "");

    for t in &targets {
        let res = if mode == "scalar" {
            viterbi_scalar(&hmm, &t.seq)
        } else {
            viterbi_striped(&hmm, &t.seq)
        };

        let evalue = compute_evalue(&hmm, res.bit_score, targets.len());
        if evalue <= max_evalue {
            println!(
                "{:20} {:>10} {:>10.1} {:>12.2e} {:>8} {:>8} {:>8} {:>8}",
                t.name, t.seq.len(), res.bit_score, evalue, res.model_start, res.model_end, res.query_start, res.query_end
            );

            if let Some(ref mut f) = domtbl_file {
                writeln!(
                    f, "{}\t{}\t{}\t{}\t{:.1}\t{:.2e}\t{}\t{}\t{}\t{}",
                    t.name, t.seq.len(), hmm.name, hmm.length, res.bit_score, evalue,
                    res.model_start, res.model_end, res.query_start, res.query_end
                ).unwrap();
            }

            if show_alignment {
                let state_str: String = res.state_path.into_iter().collect();
                println!("  States: {}", state_str);
            }
        }
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        print_usage();
        process::exit(1);
    }

    match args[1].as_str() {
        "build" => cmd_build(&args[2..]),
        "search" => cmd_search(&args[2..]),
        _ => {
            eprintln!("Unknown command: {}", args[1]);
            print_usage();
            process::exit(1);
        }
    }
}
