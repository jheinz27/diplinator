use clap::Parser;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

#[derive(Parser, Debug)]
#[command(
    name = "separate_haps_fasta",
    about = "Split a FASTA file into two haplotype files based on header tags",
    version
)]
struct Cli {
    #[arg(value_name = "INPUT", help = "input FASTA file (.fa, .fasta, .fa.gz, .fasta.gz)")]
    input: String,

    #[arg(
        short = '1',
        long,
        value_name = "TAG",
        default_value = "hap1",
        help = "tag in header identifying first haplotype"
    )]
    tag1: String,

    #[arg(
        short = '2',
        long,
        value_name = "TAG",
        default_value = "hap2",
        help = "tag in header identifying second haplotype"
    )]
    tag2: String,

    #[arg(
        short,
        long,
        value_name = "PREFIX",
        help = "output file prefix [default: derived from input filename]"
    )]
    output: Option<String>,

    #[arg(long, value_name = "FILE", help = "output file name for tag1 sequences (overrides -o/--output)")]
    out1: Option<String>,

    #[arg(long, value_name = "FILE", help = "output file name for tag2 sequences (overrides -o/--output)")]
    out2: Option<String>,

    #[arg(
        short,
        long,
        value_name = "FILE",
        help = "write sequences matching neither tag to FILE (discarded if not set)"
    )]
    unmatched: Option<String>,
}

enum Writer {
    Plain(BufWriter<File>),
    Gz(GzEncoder<BufWriter<File>>),
}

impl Write for Writer {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            Writer::Plain(w) => w.write(buf),
            Writer::Gz(w) => w.write(buf),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            Writer::Plain(w) => w.flush(),
            Writer::Gz(w) => w.flush(),
        }
    }
}

fn create_writer(path: &str, gzip: bool) -> Writer {
    let file = BufWriter::new(File::create(path).unwrap_or_else(|e| {
        eprintln!("Failed to create {}: {}", path, e);
        std::process::exit(1);
    }));
    if gzip {
        Writer::Gz(GzEncoder::new(file, Compression::default()))
    } else {
        Writer::Plain(file)
    }
}

fn main() {
    let args = Cli::parse();

    let input_path = &args.input;
    let is_gzipped = input_path.ends_with(".gz");

    // Derive the base extension (fa/fasta) and whether output should be gzipped
    let path = Path::new(input_path);
    let (stem, base_ext) = if is_gzipped {
        // e.g. genome.fa.gz -> stem=genome, ext=fa
        let without_gz = path.file_stem().unwrap().to_str().unwrap();
        let inner = Path::new(without_gz);
        let stem = inner.file_stem().unwrap_or(inner.as_os_str()).to_str().unwrap();
        let ext = inner.extension().map(|e| e.to_str().unwrap()).unwrap_or("fa");
        (stem.to_string(), ext.to_string())
    } else {
        let stem = path.file_stem().expect("No file stem").to_str().unwrap().to_string();
        let ext = path.extension().map(|e| e.to_str().unwrap()).unwrap_or("fa").to_string();
        (stem, ext)
    };

    let prefix = match &args.output {
        Some(p) => p.clone(),
        None => stem.clone(),
    };

    let gz_suffix = if is_gzipped { ".gz" } else { "" };
    let tag1_path = match &args.out1 {
        Some(p) => p.clone(),
        None => format!("{}.{}.{}{}", prefix, args.tag1, base_ext, gz_suffix),
    };
    let tag2_path = match &args.out2 {
        Some(p) => p.clone(),
        None => format!("{}.{}.{}{}", prefix, args.tag2, base_ext, gz_suffix),
    };

    // Open reader with gzip detection
    let file = File::open(input_path).expect("Failed to open input FASTA");
    let reader: Box<dyn BufRead> = if is_gzipped {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut tag1_writer = create_writer(&tag1_path, tag1_path.ends_with(".gz"));
    let mut tag2_writer = create_writer(&tag2_path, tag2_path.ends_with(".gz"));

    let mut unmatched_writer = args.unmatched.as_ref().map(|p| {
        create_writer(p, p.ends_with(".gz"))
    });

    let mut current_dest: Option<u8> = None;
    let mut tag1_count: u64 = 0;
    let mut tag2_count: u64 = 0;
    let mut unmatched_count: u64 = 0;

    for line in reader.lines() {
        let line = line.expect("Failed to read line");
        if line.starts_with('>') {
            if line.contains(args.tag1.as_str()) {
                current_dest = Some(1);
                tag1_count += 1;
            } else if line.contains(args.tag2.as_str()) {
                current_dest = Some(2);
                tag2_count += 1;
            } else {
                unmatched_count += 1;
                if unmatched_writer.is_some() {
                    current_dest = Some(0);
                } else {
                    current_dest = None;
                    eprintln!("Warning: skipping unmatched header: {}", line);
                }
            }
        }
        match current_dest {
            Some(1) => writeln!(tag1_writer, "{}", line).unwrap(),
            Some(2) => writeln!(tag2_writer, "{}", line).unwrap(),
            Some(0) => {
                if let Some(ref mut w) = unmatched_writer {
                    writeln!(w, "{}", line).unwrap();
                }
            }
            _ => {}
        }
    }

    tag1_writer.flush().unwrap();
    tag2_writer.flush().unwrap();
    if let Some(ref mut w) = unmatched_writer {
        w.flush().unwrap();
    }

    eprintln!("Done.");
    eprintln!("  {} -> {} sequence(s)", tag1_path, tag1_count);
    eprintln!("  {} -> {} sequence(s)", tag2_path, tag2_count);
    if unmatched_count > 0 {
        match &args.unmatched {
            Some(p) => eprintln!("  {} -> {} unmatched sequence(s)", p, unmatched_count),
            None => eprintln!("  {} unmatched sequence(s) discarded", unmatched_count),
        }
    }
}
