pub mod cli;
pub use cli::Cli;
pub mod paf;
pub mod sam;
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
//modeled on BWA-MEM's mem_approx_mapq_se (bwamem.c)
//confidence measure that a read is assigned to the correct haplotype
pub fn compute_hapq(score_winner: f32, score_loser: f32, n_splits: u32, match_sc: f32) -> u8 {
    if score_winner <= 0.0 {
        return 0;
    }
    //approximately the difference in matching bases btwn alignment1 and alignment2
    let diff = (score_winner - score_loser) / match_sc;
    //penalize reads with more that 3 split aligments (likely a complex region)
    let pen_split = if n_splits <= 3 { 1.0 } else { 3.0 / n_splits as f32 };
    let score = 6.02 * diff * pen_split;
    score.clamp(0.0, 60.0) as u8
    
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

//print end of run summary statistics to stderr, shared by the SAM and PAF paths
//counts and bases are both reported per category since read counts alone can mislead:
//a category can be inflated by many very short reads
//order of both arrays is [asm1, asm2, equal, unmapped]
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
    //report bases as gigabases, keeping every digit: integer division for the whole part
    //and the remainder as the 9 decimal places, so no rounding or float error is introduced
    let gbp = |n: u64| format!("{}.{:09}", n / 1_000_000_000, n % 1_000_000_000);
    //widths of the read count and gigabase columns
    let cw = total.to_string().len();
    let bw = gbp(total_bases).len();

    for (i, label) in labels.iter().take(4).enumerate() {
        eprintln!("{:<lw$} {:>cw$} reads ({:>4.1}%) ; {:>bw$} Gbps ({:>4.1}%)",
            label, counts[i], pct(counts[i], total), gbp(bases[i]), pct(bases[i], total_bases));
    }
    //pad where the percentages would be so the totals line up with the rows above
    eprintln!("{:<lw$} {:>cw$} reads          ; {:>bw$} Gbps", labels[4], total, gbp(total_bases));
}
 