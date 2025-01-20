import pandas as pd
import pybedtools
from pathlib import Path
from dotenv import load_dotenv
import os

load_dotenv()

OUTPUT_PATH = Path(os.getenv("OUTPUT_PATH"))/'garcia_ATAC/annotation'
PEAKS_N_COL = 5

def extract_tss_bed(gtf_file, tss_bed_file):
    """
    Extract TSS (start for '+' genes, end for '-' genes).
    We store: chr, tss_start, tss_end, gene_name, strand
    where tss_end = tss_start + 1
    """
    out_lines = []
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            feature_type = parts[2]
            if feature_type.lower() == "gene": 
                chrom = parts[0]
                strand = parts[6]
                start = int(parts[3]) - 1  # 0-based
                end   = int(parts[4])
                attr  = parts[8]
                
                # Extract gene name
                gene_name = ""
                if 'gene_name' in attr:
                    gene_name = attr.split('gene_name "')[1].split('"')[0]
                elif 'gene_id' in attr:
                    gene_name = attr.split('gene_id "')[1].split('"')[0]

                if strand == "+":
                    tss_pos = start  # gene start is TSS
                else:
                    tss_pos = end - 1  # gene end - 1 for 0-based

                out_lines.append(f"{chrom}\t{tss_pos}\t{tss_pos+1}\t{gene_name}\t.\t{strand}\n")

    with open(tss_bed_file, 'w') as fout:
        fout.writelines(out_lines)

    print(f"Created TSS BED: {tss_bed_file} (lines={len(out_lines)})")


def create_promoter_bed(tss_bed, promoter_bed, genome_file,
                        plus_range=(-1000, 100), minus_range=(-100, 1000)):
    """
    Create promoter regions from the TSS bed. For + strand: (tss-1000, tss+100).
    For - strand: (tss-100, tss+1000). We use pybedtools.slop().
    """
    tss = pybedtools.BedTool(tss_bed)
    # Separate by strand
    plus  = tss.filter(lambda x: x.strand == "+").saveas("tmp.tss.plus.bed")
    minus = tss.filter(lambda x: x.strand == "-").saveas("tmp.tss.minus.bed")

    # plus strand: left=1000, right=100
    plus_slop = plus.slop(
        g=genome_file,
        l=abs(plus_range[0]),
        r=plus_range[1],
        s=True
    )

    # minus strand: left=100, right=1000
    minus_slop = minus.slop(
        g=genome_file,
        l=abs(minus_range[0]),
        r=minus_range[1],
        s=True
    )

    # Combine
    merged = plus_slop.cat(minus_slop, postmerge=False).saveas(promoter_bed)
    print(f"Created promoter BED: {promoter_bed}, lines={len(merged)}")

    # cleanup
    os.remove("tmp.tss.plus.bed")
    os.remove("tmp.tss.minus.bed")


def create_gene_body_bed(gtf_file, body_bed_file):
    """
    Create a BED representing the entire gene body (start..end).
    We'll store: chr, start, end, gene_name, strand
    """
    out_lines = []
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            feature_type = parts[2]
            if feature_type.lower() == "gene":
                chrom = parts[0]
                strand = parts[6]
                start = int(parts[3]) - 1  # 0-based
                end   = int(parts[4])
                attr  = parts[8]
                
                # Extract gene name
                gene_name = ""
                if 'gene_name' in attr:
                    gene_name = attr.split('gene_name "')[1].split('"')[0]
                elif 'gene_id' in attr:
                    gene_name = attr.split('gene_id "')[1].split('"')[0]

                out_lines.append(f"{chrom}\t{start}\t{end}\t{gene_name}\t.\t{strand}\n")

    with open(body_bed_file, 'w') as fout:
        fout.writelines(out_lines)

    print(f"Created gene-body BED: {body_bed_file}, lines={len(out_lines)}")


def subtract_bed(a_bed, b_bed, out_bed):
    """
    'subtract' (bedtools intersect -v). Keep intervals in a_bed that do NOT intersect b_bed.
    """
    a = pybedtools.BedTool(a_bed)
    b = pybedtools.BedTool(b_bed)
    result = a.intersect(b, v=True)
    result.saveas(out_bed)
    print(f"{len(result)} remaining. ")

    return out_bed


def annotate_promoter_peaks(peaks_bed, promoter_bed, out_bed):
    """
    Return intersections of peaks with promoters. 
    """
    peaks = pybedtools.BedTool(peaks_bed)
    print(f"Total number of peaks to annotate: {len(peaks)}")
    promos = pybedtools.BedTool(promoter_bed)
    inter = peaks.intersect(promos, wa=True, wb=True)
    inter.saveas(out_bed)
    print(f"Created annotated promoter peaks BED: {out_bed}, lines={len(inter)}")
    return out_bed


def annotate_distal_peaks(peaks_bed, tss_bed, out_bed, cutoff=200000):
    """
    Find peaks within ±200 kb of the nearest TSS (bedtools closest -D ref).
    Output only those that meet the cutoff.
    """
    peaks = pybedtools.BedTool(peaks_bed)
    tss   = pybedtools.BedTool(tss_bed)
    # bedtools closest, last column = distance
    close = peaks.sort().closest(tss.sort(), D="ref")
    
    # Filter by distance
    def keep_distal(feature):
        dist = int(feature[-1])
        if abs(dist) <= cutoff:
            return True
        return False
    
    close.filter(keep_distal).saveas(out_bed)
    return out_bed


def annotate_body_peaks(peaks_bed, body_bed, out_bed):
    """
    Peaks that overlap gene bodies (bedtools intersect wa wb).
    """
    peaks = pybedtools.BedTool(peaks_bed)
    bodies = pybedtools.BedTool(body_bed)
    inter = peaks.intersect(bodies, wa=True, wb=True)
    inter.saveas(out_bed)
    return out_bed


def label_peaks_with_category(
    promoter_intersections, distal_intersections, body_intersections, intergenic_intersections,
    all_peaks_bed, final_csv):
    """
    Combine all categories into a single CSV with columns:
       peak_chr, peak_start, peak_end, gene_name, annotation, distance?
    We will do a simple approach:
       - If peak found in promoter_intersections => 'promoter'
       - Else if peak found in distal_intersections => 'distal'
       - Else if peak found in body_intersections => 'distal_body'
       - Else => 'intergenic'
       
    """
    # Read all peaks as a DF
    all_peaks_df = pd.read_csv(all_peaks_bed, sep="\t", header=None)
    all_peaks_df.columns = ["peak_chr","peak_start","peak_end","peak_name","peak_score"]
    
    # For each category, gather the peak IDs
    def read_intersections(inter_bed):
        df = pd.read_csv(inter_bed, sep="\t", header=None)
         
        df.columns = ["peak_chr","peak_start","peak_end","peak_name","peak_score", 
                        "gene_chr","gene_promoter_start","gene_promoter_end","gene_name", "gene_score", "gene_strand", "distance"][:len(df.columns)]                
        return df
    
    df_promoter = read_intersections(promoter_intersections)
    df_promoter['annotation_type'] = "promoter_peak"
    df_distal   = read_intersections(distal_intersections)
    df_distal['annotation_type'] = "distal_peak"
    df_body     = read_intersections(body_intersections)
    df_body['annotation_type'] = "distal_peak_body"    
    df_intergenic = read_intersections(intergenic_intersections)
    df_intergenic['annotation_type'] = "intergenic"    
    
    final_df = pd.concat([df_promoter, df_distal, df_body, df_intergenic])
    print(final_df.annotation_type.value_counts())
    final_df.to_csv(final_csv, index=False)
    print(f"Annotation categories saved to {final_csv}")
    # Sanity checks
    peaks_unique = final_df.iloc[:, :PEAKS_N_COL].drop_duplicates()
    print(f'Matched {len(peaks_unique)} vs {len(all_peaks_df)}')
    sorted_peaks_unique = peaks_unique.sort_values(by=list(peaks_unique.columns)).reset_index(drop=True)
    sorted_all_peaks_df = all_peaks_df.sort_values(by=list(all_peaks_df.columns)).reset_index(drop=True)

    if not sorted_peaks_unique.equals(sorted_all_peaks_df):
        raise ValueError("peaks_unique and all_peaks_df are different (ignoring headers and order).")
    else:
        print("peaks_unique and all_peaks_df are the same (ignoring headers and order).")

###############################################################################
# Main advanced annotation logic
###############################################################################
def run_strand_aware_annotation(
    gtf_file,
    peaks_bed,
    genome_file,
    output_prefix,
    promoter_plus_minus = ((-1000, 100), (-100, 1000)),
    distal_cutoff=200000
):
    """
    Full pipeline:
      1) Extract TSS
      2) Create promoter bed
      3) Annotate promoter peaks
      4) Subtract promoter from peaks => "non-promoter"
      5) Distal annotation (closest TSS within ±200 kb)
      6) Subtract distal => "non-promoter-non-distal"
      7) Gene body intersection
      8) Subtract body => intergenic
      9) Combine final labels in a CSV
    """
    print("==> 1) Extract TSS from GTF...")
    tss_bed = f"{output_prefix}.tss.bed"
    extract_tss_bed(gtf_file, tss_bed)

    print("==> 2) Create promoter bed (strand-aware)...")
    promoter_bed = f"{output_prefix}.promoters.bed"
    create_promoter_bed(
        tss_bed, promoter_bed, genome_file,
        plus_range=promoter_plus_minus[0],
        minus_range=promoter_plus_minus[1],
    )

    print("==> 3) Annotate promoter peaks...")
    promoter_hits_bed = f"{output_prefix}.promoter_hits.bed"
    annotate_promoter_peaks(peaks_bed, promoter_bed, promoter_hits_bed)

    # Unique promoter peaks
    unique_promoter_peaks = f"{output_prefix}.promoter_peaks.unique.bed"
    # We only need the peak portion from the intersection (first columns)
    promoter_intersect_df = pd.read_csv(promoter_hits_bed, sep="\t", header=None)
  
    peak_cols = list(range(PEAKS_N_COL))
    promoter_peaks_unique = promoter_intersect_df.iloc[:, peak_cols].drop_duplicates()
    promoter_peaks_unique.to_csv(unique_promoter_peaks, sep="\t", header=False, index=False)
    print(f"{len(promoter_peaks_unique)} were matched using the 1st overlap rule. ")

    print("==> 4) Remove promoter peaks => non-promoter peaks")
    non_promoter_bed = f"{output_prefix}.nonpromoter.bed"
    subtract_bed(peaks_bed, unique_promoter_peaks, non_promoter_bed)

    print("==> 5) Distal annotation (±200 kb from TSS)")
    distal_bed = f"{output_prefix}.distal.bed"
    annotate_distal_peaks(non_promoter_bed, tss_bed, distal_bed, cutoff=distal_cutoff)

    # unique distal peaks
    distal_df = pd.read_csv(distal_bed, sep="\t", header=None)
    distal_peaks = distal_df.iloc[:, :PEAKS_N_COL].drop_duplicates()
    distal_peaks_bed = f"{output_prefix}.distal_peaks.unique.bed"
    distal_peaks.to_csv(distal_peaks_bed, sep="\t", header=False, index=False)

    # leftover after distal
    non_distal_bed = f"{output_prefix}.nonpromoter_nondistal.bed"
    subtract_bed(non_promoter_bed, distal_peaks_bed, non_distal_bed)

    print("==> 6) Gene body annotation (for leftover peaks)")
    body_bed = f"{output_prefix}.gene_body.bed"
    create_gene_body_bed(gtf_file, body_bed)

    body_hits_bed = f"{output_prefix}.body_hits.bed"
    annotate_body_peaks(non_distal_bed, body_bed, body_hits_bed)

    # leftover => intergenic
    intergenic_bed = f"{output_prefix}.intergenic.bed"
    subtract_bed(non_distal_bed, body_hits_bed, intergenic_bed)

    print("==> 7) Combine annotation categories into a single CSV")
    final_csv = f"{output_prefix}.final_annotation.csv"
    label_peaks_with_category(
        promoter_intersections=promoter_hits_bed,
        distal_intersections=distal_bed,
        body_intersections=body_hits_bed,
        intergenic_intersections=intergenic_bed,
        all_peaks_bed=peaks_bed,
        final_csv=final_csv
    )
    print("Done.")


###############################################################################
# Main script
###############################################################################
if __name__ == "__main__":
    # Inputs
    gtf_file = "/home/bdobre/resources/gencode.v47.basic.annotation.gtf"
    peaks_file = "/mnt/storage/outputs/garcia_ATAC/atac_preprocessing_combined/consensus_peak_calling/consensus_regions.bed"
    genome_file = "/home/bdobre/resources/hg38.genome"  # A file with chrom sizes: chr <tab> length

    # Output prefix for all intermediate files
    output_prefix = str(OUTPUT_PATH / "strand_aware_annotation")
    
    # Run advanced logic
    run_strand_aware_annotation(
        gtf_file=gtf_file,
        peaks_bed=peaks_file,
        genome_file=genome_file,
        output_prefix=output_prefix,
        promoter_plus_minus=((-1000, 100), (-100, 1000)),  # (plus_range, minus_range)
        distal_cutoff=200000
    )
