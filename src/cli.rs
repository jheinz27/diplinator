use clap::{Parser, ValueEnum};


#[derive(Parser, Debug)]
#[command( name = "HipHap", about = "HipHap: Choose the best alignment to each haploid of a diploid assembly", version)]
pub struct Cli {
    //
    #[arg(value_name = "ASM1", help="asm1 alignment file (sam/bam/cram/paf)")]
    pub asm1: String,

    #[arg(value_name = "ASM2", help="asm2 alignment file (sam/bam/cram/paf)")]
    pub asm2: String,

    #[arg(short='1', long, value_name = "NAME", default_value = "asm1", help="label for asm1 sample (used in output file names and summary)")]
    pub s1: String,

    #[arg(short='2', long, value_name = "NAME", default_value = "asm2", help="label for asm2 sample (used in output file names and summary)")]
    pub s2: String,

    // inputs are PAF files
    #[arg(long, default_value_t = false, help = "input files are PAF")]
    pub paf: bool,

    //use ms score rather than AS score
    #[arg(long, default_value_t = false, help = "use ms:i: tag rather than AS:i: for alignment score")]
    pub ms: bool,

    // write tied reads to both output files (requires -p: two primaries cannot share one file)
    #[arg(short, long, default_value_t = false, help = "write reads with equal alignment scores to both output files (requires -p)")]
    pub both: bool,

    // write one file per haplotype rather than a single merged output file (merged is the default)
    #[arg(short = 'p', long, default_value_t = false, help = "write one file per haplotype instead of a single merged output file")]
    pub partition: bool,

    // output path: the merged file name, or the shared stem for the pair under --partition
    #[arg(short = 'o', long, value_name = "FILE", help = "output file [default: hiphap_{s1}_{s2}_merged.*]; with -p this is the stem for {out}_{s1}.* and {out}_{s2}.*")]
    pub output: Option<String>,

    // combined reference FASTA for writing a merged CRAM (must contain all contigs of both haplotypes)
    #[arg(long, value_name = "FILE", required = false, help = "combined reference FASTA for merged CRAM output (must contain all contigs of both inputs); required for merged CRAM output")]
    pub ref_merged: Option<String>,

    // where to write reads unmapped in both assemblies
    #[arg(short, long, value_name = "DEST", default_value = "asm1", help="where to write reads unmapped in both assemblies: asm1, asm2, or discard")]
    pub unmapped: UnmappedDest,

    #[arg(long, value_name = "FILE", required = false, help="reference FASTA for cram file (asm1)")]
    pub ref1: Option<String>,

    #[arg(long, value_name = "FILE", required = false, help="reference FASTA for cram file (asm2)")]
    pub ref2: Option<String>,

    // per-base match score from aligner scoring scheme (used in HAPQ calculation)
    #[arg(short = 'A' , long, value_name = "FLOAT", help = "per-base match score from aligner scoring scheme (auto-estimated from ms:i tags if omitted)")]
    pub match_sc: Option<f32>,

    // skip HAPQ score calculation and hq tag output (for non-haplotype comparisons)
    #[arg(long, default_value_t = false, help = "skip HAPQ score calculation and hq tag output (e.g. for comparing GRCh38 vs CHM13)")]
    pub no_hapq: bool,

    // disable writing the list of chromosome-spanning reads
    #[arg(long, default_value_t = false, help = "disable writing the chromosome-spanning reads file (*_span_chrom.fastq, or .txt for PAF)")]
    pub no_span_chrom: bool,

    // number of total threads to use; the default follows the mode, so resolve it with
    // Cli::resolved_threads() rather than reading this field directly
    #[arg(short, long,value_name = "INT", help = "total thread pool size [default: 6; 8 with -p]. For BAM/CRAM output the writer gets 4x a reader when merging and 3x with -p; 1x for SAM")]
    pub threads: Option<usize>
}

impl Cli {
    //the thread budget, defaulting to the smallest one that fits the mode at the intended
    //write:read ratio: merged runs 2 readers + 1 writer (2 + 4 = 6), -p runs 2 readers +
    //2 writers (2 + 3 + 3 = 8). See plan_threads for how the budget is then divided.
    pub fn resolved_threads(&self) -> usize {
        self.threads.unwrap_or(if self.partition { 8 } else { 6 })
    }
}

#[derive(Debug, Clone, ValueEnum)]
pub enum UnmappedDest {
    Asm1,
    Asm2,
    Discard,
}
