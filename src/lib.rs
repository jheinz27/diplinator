pub mod cli;
pub use cli::Cli;
pub mod paf;
pub mod sam;
use std::cmp::max;
use std::hash::{Hash, Hasher};
use twox_hash::XxHash64;

//enum to store best alignment of read
pub enum Winner {
    Asm1,
    Asm2,
    Both,
    Unmapped,
}

//if read has identical alignment to both haps,
//chose which hap to report randomly with equal likelihoods
//use last bit of hash of read ID (as bytes) as random assignment
pub fn choose_random(id: &[u8]) -> Winner {
    //XxHash64 provides reproducible assignment bc is deterministic
    let mut hasher = XxHash64::with_seed(42);
    id.hash(&mut hasher);
    if hasher.finish() & 1 == 0 { Winner::Asm1 } else { Winner::Asm2 }
}

//compute haplotype assignment quality (HAPQ) score
//confidence measure that a read is assigned to the correct haplotype
pub fn compute_hapq(score_winner: f32, score_loser: f32, n_splits: u32, match_sc: f32) -> u8 {
    
    //scores are identical between haploypes: hapq = 0 
    if score_winner <= score_loser {
        return 0;
    }
    //approximately the difference in matching bases btwn alignment1 and alignment2
    let diff = (score_winner - score_loser) / match_sc;
    //penalize reads with more that 3 split aligments (likely a complex region)
    let pen_split = if n_splits <= 3 { 1.0 } else { 3.0 / n_splits as f32 };
    let score = 6.02 * diff * pen_split;
    //hapq = 0 means that the AS scores were truly identical, special case 
    //a very low hapq score (<1) arising from the emperical penalty gets hapq=1
    (score.clamp(0.0, 60.0) as u8).max(1)
    
}

//helper function to merge any read alignment segments that overlap in read coordinates
//returns count of unique bps of the read contained in any alignment segment
pub fn merge_intervals(intervals: &mut [(u32, u32)]) -> u32 {
    //sort cluster by read start location of alignment segment
    intervals.sort_unstable_by_key(|k| k.0);
    
    let mut read_bps_aligned = 0; 
    if !intervals.is_empty() {
        //initialize at first interval
        let (mut cur_start, mut cur_end) = intervals[0]; 
        //iterate through intervals and merge adjacent overlapping intervals
        for &(next_start, next_end) in intervals.iter().skip(1) { 
            if next_start < cur_end {
                // intervals overlap, so extend
                if next_end > cur_end {
                    cur_end = next_end;
                }
            } else {
                //no further overlap, add length and start over with next interval grouping
                read_bps_aligned += cur_end - cur_start;
                cur_start = next_start;
                cur_end = next_end;
            }
        
        }
        //add final overlap segment
        read_bps_aligned += cur_end - cur_start
    }
    read_bps_aligned
}

//strip a trailing alignment-file extension so -o can seed derived file names
//(the span-chrom file), and so the extension can be compared against the output format
pub fn strip_aln_ext(path: &str) -> &str {
    for e in [".sam", ".bam", ".cram", ".paf"] {
        if let Some(stem) = path.strip_suffix(e) {
            return stem;
        }
    }
    path
}

//resolve every output path for both modes, shared by the SAM and PAF backends so the two
//cannot drift apart. ext is the alignment extension taken from the input format
//(".sam"/".bam"/".cram"/".paf"); span_ext is ".fastq" for SAM/BAM/CRAM and ".txt" for PAF.
//-> (primary, secondary, span_chrom); secondary is None in merged mode
pub fn output_paths(args: &Cli, ext: &str, span_ext: &str) -> (String, Option<String>, String) {
    //the span file follows the alignment output, so with -o it lands in the same place
    //rather than in the working directory under an unrelated name
    let span_stem = match &args.output {
        Some(o) => strip_aln_ext(o).to_string(),
        None => format!("hiphap_{}_{}", args.s1, args.s2),
    };
    let span = format!("{}_span_chrom{}", span_stem, span_ext);

    if args.partition {
        //-o is the stem the two per-haplotype files share; without it the historical
        //hiphap_{s1}{ext} / hiphap_{s2}{ext} names are unchanged. Strip any extension the
        //user supplied so '-o sample.bam' gives sample_mat.bam, not sample.bam_mat.bam
        let stem = match &args.output {
            Some(o) => strip_aln_ext(o).to_string(),
            None => "hiphap".to_string(),
        };
        (
            format!("{}_{}{}", stem, args.s1, ext),
            Some(format!("{}_{}{}", stem, args.s2, ext)),
            span,
        )
    } else {
        //-o names the merged file outright
        let merged = args
            .output
            .clone()
            .unwrap_or_else(|| format!("hiphap_{}_{}_merged{}", args.s1, args.s2, ext));
        (merged, None, span)
    }
}

//the output format always follows the input format, so an -o extension that disagrees is
//silently ignored; say so rather than leaving SAM text in a file named .bam
pub fn warn_output_ext_mismatch(args: &Cli, ext: &str) {
    let Some(o) = &args.output else { return };
    let stem = strip_aln_ext(o);
    //no recognised extension to disagree with
    if stem.len() == o.len() {
        return;
    }
    let given = &o[stem.len()..];
    if !given.eq_ignore_ascii_case(ext) {
        eprintln!(
            "Warning: output format is {} (taken from the input); the '{}' extension of '{}' is ignored",
            ext.trim_start_matches('.').to_uppercase(), given, o
        );
    }
}

//how a thread budget is divided between the htslib readers and writers
pub struct ThreadPlan {
    //threads for each reader (there are always two, one per assembly)
    pub reader: usize,
    //threads for each writer (the single merged writer, or one per haplotype with -p)
    pub writer: usize,
}

//divide the --threads budget between the readers and the writers. A reader counts as one share
//and a writer as writer_weight shares: for BGZF/CRAM output 4 merged and 3 per writer with -p,
//where compression costs several times the matching decompression, and 1 for uncompressed SAM
//text, where there is nothing to compress. With two readers that makes the natural budget
//2 + 4 = 6 merged and 2 + 3 + 3 = 8 with -p. A request too small to give every file one thread
//is raised silently, and a leftover thread that will not divide evenly between two writers is
//left idle so the pair stays symmetric.
pub fn plan_threads(
    requested: usize,
    n_readers: usize,
    n_writers: usize,
    writer_weight: usize,
) -> ThreadPlan {
    debug_assert!(n_readers > 0 && n_writers > 0 && writer_weight > 0);

    //every reader and every writer needs a thread of its own
    let total = max(requested, n_readers + n_writers);
    let unit = n_readers + n_writers * writer_weight;

    //the reader share of the budget, rounded half up so readers grow smoothly with the budget
    //instead of only stepping at exact multiples of unit
    let mut reader = (2 * total + unit) / (2 * unit);
    //but never so much that a writer is left with nothing
    reader = reader.min((total - n_writers) / n_readers).max(1);
    let mut writer = (total - n_readers * reader) / n_writers;

    //for compressed output the writer is the expensive side, so it must never end up with fewer
    //threads than a reader (only reachable at tiny budgets). Plain SAM is left alone: its target
    //is 1:1, and e.g. 5 threads merged fit 2/2/1 more closely than 1/1/3.
    while writer_weight > 1 && writer < reader && reader > 1 {
        reader -= 1;
        writer = (total - n_readers * reader) / n_writers;
    }

    ThreadPlan { reader, writer }
}

//print end of run summary statistics to stderr, shared by the SAM and PAF paths

pub fn print_summary(s1: &str, s2: &str, counts: [u64; 4], bases: [u64; 4]) {
    let total: u64 = counts.iter().sum();
    let total_bases: u64 = bases.iter().sum();
    //avoid NaN% when nothing was parsed (e.g. empty inputs)
    let pct = |n: u64, denom: u64| if denom == 0 { 0.0 } else { n as f64 / denom as f64 * 100.0 };

    //build labels up front so the number columns line up for any sample name length
    let labels = [
        format!("Reads aligned better to {}:", s1),
        format!("Reads aligned better to {}:", s2),
        "Reads with equal scores:".to_string(),
        "Reads unmapped to both:".to_string(),
        "Total reads parsed:".to_string(),
    ];
    let lw = labels.iter().map(|l| l.len()).max().unwrap_or(0);
    //report bases as gigabases
    let gbp = |n: u64| format!("{}.{:09}", n / 1_000_000_000, n % 1_000_000_000);
    let cw = total.to_string().len();
    let bw = gbp(total_bases).len();

    for (i, label) in labels.iter().take(4).enumerate() {
        eprintln!("{:<lw$} {:>cw$} reads ({:>4.1}%) ; {:>bw$} Gbps ({:>4.1}%)",
            label, counts[i], pct(counts[i], total), gbp(bases[i]), pct(bases[i], total_bases));
    }
    //pad where the percentages would be so the totals line up with the rows above
    eprintln!("{:<lw$} {:>cw$} reads          ; {:>bw$} Gbps", labels[4], total, gbp(total_bases));
}

#[cfg(test)]
mod tests {
    use super::plan_threads;

    //the four shapes hiphap actually runs: two readers, one or two writers, compressed or not
    const MERGED_COMPRESSED: (usize, usize) = (1, 4);
    const MERGED_SAM: (usize, usize) = (1, 1);
    const PART_COMPRESSED: (usize, usize) = (2, 3);
    const PART_SAM: (usize, usize) = (2, 1);

    //(reader, writer) for a given budget and shape
    fn split(requested: usize, shape: (usize, usize)) -> (usize, usize) {
        let p = plan_threads(requested, 2, shape.0, shape.1);
        (p.reader, p.writer)
    }

    #[test]
    fn matches_the_default_budgets() {
        //6 merged: one thread per reader, four for the single writer
        assert_eq!(split(6, MERGED_COMPRESSED), (1, 4));
        //8 partitioned: one thread per reader, three for each of the two writers
        assert_eq!(split(8, PART_COMPRESSED), (1, 3));
    }

    #[test]
    fn raises_budgets_below_one_thread_per_file() {
        for t in 0..=3 {
            assert_eq!(split(t, MERGED_COMPRESSED), (1, 1), "merged -t {}", t);
        }
        for t in 0..=4 {
            assert_eq!(split(t, PART_COMPRESSED), (1, 1), "partition -t {}", t);
        }
    }

    #[test]
    fn merged_compressed_table() {
        let expected = [
            (3, (1, 1)), (4, (1, 2)), (5, (1, 3)), (6, (1, 4)),
            (8, (1, 6)), (10, (2, 6)), (16, (3, 10)), (32, (5, 22)),
        ];
        for (t, want) in expected {
            assert_eq!(split(t, MERGED_COMPRESSED), want, "-t {}", t);
            //a lone writer takes whatever the readers leave, so nothing ever idles
            assert_eq!(2 * want.0 + want.1, t, "-t {} leaves a thread idle", t);
        }
    }

    #[test]
    fn partition_compressed_table() {
        //(budget, (reader, writer), threads left idle)
        let expected = [
            (4, (1, 1), 0), (5, (1, 1), 1), (6, (1, 2), 0), (7, (1, 2), 1),
            (8, (1, 3), 0), (9, (1, 3), 1), (12, (2, 4), 0), (16, (2, 6), 0),
            (20, (3, 7), 0), (24, (3, 9), 0), (32, (4, 12), 0),
        ];
        for (t, want, idle) in expected {
            assert_eq!(split(t, PART_COMPRESSED), want, "-t {}", t);
            assert_eq!(t - (2 * want.0 + 2 * want.1), idle, "-t {} idle count", t);
        }
    }

    #[test]
    fn sam_splits_evenly() {
        //exact 1:1 whenever the budget divides by the number of files
        assert_eq!(split(3, MERGED_SAM), (1, 1));
        assert_eq!(split(6, MERGED_SAM), (2, 2));
        assert_eq!(split(9, MERGED_SAM), (3, 3));
        assert_eq!(split(8, PART_SAM), (2, 2));
        assert_eq!(split(12, PART_SAM), (3, 3));
        //and the closest fit otherwise, ties going to the readers
        assert_eq!(split(5, MERGED_SAM), (2, 1));
        assert_eq!(split(6, PART_SAM), (2, 1));
    }

    #[test]
    fn invariants_hold_across_every_budget() {
        for shape in [MERGED_COMPRESSED, MERGED_SAM, PART_COMPRESSED, PART_SAM] {
            let (n_writers, weight) = shape;
            for requested in 0..=128 {
                let (r, w) = split(requested, shape);
                //nothing is ever left without a thread
                assert!(r >= 1 && w >= 1, "{:?} -t {} gave {}/{}", shape, requested, r, w);
                //and the budget is never overspent (a too-small request is raised to the minimum)
                let budget = requested.max(2 + n_writers);
                let assigned = 2 * r + n_writers * w;
                assert!(assigned <= budget, "{:?} -t {} spent {} of {}", shape, requested, assigned, budget);
                //at most one thread idles, and only when two writers cannot split the remainder
                assert!(budget - assigned <= n_writers - 1, "{:?} -t {} idled {}", shape, requested, budget - assigned);
                //compressed output must never starve the writer relative to a reader
                if weight > 1 {
                    assert!(w >= r, "{:?} -t {} gave writer {} < reader {}", shape, requested, w, r);
                }
            }
        }
    }
}
 