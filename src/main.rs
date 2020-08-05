use std::env;
use rust_htslib::{bam, bam::Read};
use rust_htslib::bam::record::CigarStringView;
use rust_htslib::bam::record::Aux;

// We need to parse the CIGAR string to find the portion
// of the read that is actually aligned so we can save the correct
// substring of BD and BI.
fn aligned_section(c: &CigarStringView, seqlen: usize) -> (usize, usize) {
    let offset_pre = c.leading_hardclips();
    let offset_post = c.trailing_hardclips();
    return (offset_pre as usize, seqlen - (offset_post as usize));
}

fn main() {
    // first, handle the arguments
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        println!("Usage: doppelganger [input bam] > output.bam");
        return;
    };

    let path = &args[1];
    let mut bam = bam::Reader::from_path(path).unwrap();
    let header = bam::Header::from_template(bam.header());
    let mut out = bam::Writer::from_stdout(&header, bam::Format::BAM).unwrap();
    let mut record_counter = 0;

    for rec in bam.records() {
        let record = rec.unwrap();
        let mut record_out = record.clone();
        if record.is_supplementary() {
            // extract the horrible non-standard tags that sentieon
            // happily polutes the consensus bam with.
            // these causes problems when the aligned sequence is shorter
            // than the read itself, as these tags (BI and BD) are copied
            // from the bam file that was used to generate the consensus reads.
            // @see https://support.sentieon.com/appnotes/umi/
            // these tags (according to sentieon) improves variant calling quality,
            // but GATK-based tools have moved away from them (at least for data
            // generated on illumina machines), and GATK4 Mutect2 fails when
            // len(BD) != len(QUAL).
            let bi = record.aux(b"BI").unwrap().string();
            let bd = record.aux(b"BD").unwrap().string();
            let offsets = aligned_section(&record.cigar().clone(), bi.len());
            let clipped_bi = Aux::String(&bi[offsets.0..offsets.1]);
            let clipped_bd = Aux::String(&bd[offsets.0..offsets.1]);

            if !record_out.remove_aux(b"BI") || !record_out.remove_aux(b"BD") {
                panic!("doppelganger: unable to delete old BI and BD tags");
            }

            record_out.push_aux(b"BI", &clipped_bi);
            record_out.push_aux(b"BD", &clipped_bd);
            record_counter += 1;
        }
        out.write(&record_out).unwrap();
    }
    eprintln!("doppelganger processed a total of {} supplementary reads", record_counter);
}
