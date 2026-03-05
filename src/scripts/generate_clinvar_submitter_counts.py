"""
Generate a TSV file of ClinVar pathogenic submitter counts per variant.

Uses ClinVar's submission_summary.txt and variant_summary.txt (from NCBI FTP)
to count the number of independent submissions classifying each variant as
pathogenic or likely pathogenic. This provides per-submission granularity
equivalent to Nirvana's RCV-level ClinVar data.

Download the required files:
    wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/submission_summary.txt.gz
    wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz

Output format (TSV, no header):
    chrom\tpos\tref\talt\tpathogenic_submitter_count

Usage:
    python3 src/scripts/generate_clinvar_submitter_counts.py <variant_summary.txt.gz> <submission_summary.txt.gz> <output.tsv> [--assembly GRCh37]
"""
import gzip
import sys


def load_variation_id_to_coords(variant_summary_path, assembly="GRCh37"):
    """
    Load VariationID -> (chrom, pos, ref, alt) mapping from variant_summary.txt.

    Uses PositionVCF, ReferenceAlleleVCF, AlternateAlleleVCF columns for
    VCF-compatible coordinates.
    """
    variation_to_coords = {}
    open_func = gzip.open if variant_summary_path.endswith(".gz") else open

    with open_func(variant_summary_path, "rt") as f:
        header = f.readline().strip().split("\t")
        col = {name: i for i, name in enumerate(header)}

        # Identify required columns
        required = ["VariationID", "Assembly", "Chromosome", "PositionVCF",
                     "ReferenceAlleleVCF", "AlternateAlleleVCF"]
        for r in required:
            if r not in col:
                print(f"ERROR: Column '{r}' not found in variant_summary.txt header")
                print(f"Available columns: {header}")
                sys.exit(1)

        for line in f:
            fields = line.strip().split("\t")
            if len(fields) <= max(col[r] for r in required):
                continue

            if fields[col["Assembly"]] != assembly:
                continue

            var_id = fields[col["VariationID"]]
            chrom = fields[col["Chromosome"]]
            pos = fields[col["PositionVCF"]]
            ref = fields[col["ReferenceAlleleVCF"]]
            alt = fields[col["AlternateAlleleVCF"]]

            if not pos or pos == "-1" or pos == "na":
                continue

            # Add chr prefix to match VEP output
            if not chrom.startswith("chr"):
                chrom = f"chr{chrom}"

            variation_to_coords[var_id] = (chrom, pos, ref, alt)

    print(f"Loaded {len(variation_to_coords)} variant coordinates from variant_summary.txt ({assembly})")
    return variation_to_coords


def count_pathogenic_submissions(submission_summary_path):
    """
    Count pathogenic/likely pathogenic submissions per VariationID from submission_summary.txt.

    Each row is one SCV (submission). We count rows where ClinicalSignificance
    contains 'pathogenic' (covers Pathogenic, Likely pathogenic, Pathogenic/Likely pathogenic).
    """
    variation_to_patho_count = {}
    open_func = gzip.open if submission_summary_path.endswith(".gz") else open

    with open_func(submission_summary_path, "rt") as f:
        header_found = False
        var_id_idx = None
        sig_idx = None

        for line in f:
            if not header_found:
                if line.startswith("#"):
                    # Keep parsing — the last '#' line is the header
                    header = line.lstrip("#").strip().split("\t")
                    col = {name: i for i, name in enumerate(header)}
                    if "VariationID" in col and "ClinicalSignificance" in col:
                        var_id_idx = col["VariationID"]
                        sig_idx = col["ClinicalSignificance"]
                        header_found = True
                    continue
                else:
                    if var_id_idx is None:
                        print(f"ERROR: Header with VariationID and ClinicalSignificance not found")
                        sys.exit(1)
                    header_found = True

            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) <= max(var_id_idx, sig_idx):
                continue

            var_id = fields[var_id_idx]
            sig = fields[sig_idx].lower()

            if "pathogenic" in sig and "conflict" not in sig:
                variation_to_patho_count[var_id] = variation_to_patho_count.get(var_id, 0) + 1

    print(f"Found {len(variation_to_patho_count)} variants with pathogenic submissions")
    return variation_to_patho_count


def generate_submitter_counts(variant_summary_path, submission_summary_path, output_path, assembly="GRCh37"):
    """Join coordinates with submission counts and write output."""
    variation_to_coords = load_variation_id_to_coords(variant_summary_path, assembly)
    variation_to_patho_count = count_pathogenic_submissions(submission_summary_path)

    variants_written = 0
    with open(output_path, "w") as out:
        for var_id, count in variation_to_patho_count.items():
            coords = variation_to_coords.get(var_id)
            if not coords:
                continue
            chrom, pos, ref, alt = coords
            out.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{count}\n")
            variants_written += 1

    print(f"Wrote {variants_written} variants with pathogenic submitters to {output_path}")


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print(f"Usage: {sys.argv[0]} <variant_summary.txt.gz> <submission_summary.txt.gz> <output.tsv> [--assembly GRCh37]")
        sys.exit(1)

    variant_summary = sys.argv[1]
    submission_summary = sys.argv[2]
    output = sys.argv[3]
    assembly = "GRCh37"
    if "--assembly" in sys.argv:
        idx = sys.argv.index("--assembly")
        if idx + 1 < len(sys.argv):
            assembly = sys.argv[idx + 1]

    generate_submitter_counts(variant_summary, submission_summary, output, assembly)
