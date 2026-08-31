"""
Standalone gene-level data loader for BIAS-2015.

Loads gene-disease validity (ClinGen) and gene constraint (gnomAD LOEUF/pLI)
data from external files, producing an hgnc_to_gene_data dict compatible with
both the Nirvana and VEP classification pipelines.

This provides an annotator-independent source of truth for gene-level data,
replacing the annotator-specific gene sections (e.g., Nirvana's embedded genes).

Data sources:
    - ClinGen Gene-Disease Validity: https://search.clinicalgenome.org/kb/gene-validity/download
    - gnomAD Gene Constraint (v2.1.1): https://gnomad.broadinstitute.org/downloads#v2-constraint
"""
import csv
import logging


# ClinGen MOI abbreviation to full inheritance string (matching Nirvana format)
CLINGEN_MOI_MAP = {
    "AD": "Autosomal dominant",
    "AR": "Autosomal recessive",
    "XL": "X-linked",
    "SD": "Semi-dominant",
    "MT": "Mitochondrial",
    "UD": "Undetermined",
}


def load_clingen_gene_validity(clingen_csv_path):
    """
    Load ClinGen Gene-Disease Validity curations from CSV.

    Returns a dict mapping gene symbol to a list of validity entries,
    matching the Nirvana clingenGeneValidity format:
        [{diseaseId, disease, inheritance, classification, classificationDate}]

    Args:
        clingen_csv_path (str): Path to ClinGen Gene-Disease Validity CSV download

    Returns:
        dict: gene_symbol -> list of validity entries
    """
    gene_to_validity = {}
    with open(clingen_csv_path, 'r') as f:
        reader = csv.reader(f)
        for i, row in enumerate(reader):
            # Skip the 4 header/separator lines and the column header + separator lines
            if i < 6:
                continue
            if len(row) < 9:
                continue
            gene_symbol = row[0].strip()
            disease_label = row[2].strip()
            disease_id = row[3].strip()
            moi = row[4].strip()
            classification = row[6].strip().lower()
            classification_date = row[8].strip()

            inheritance = CLINGEN_MOI_MAP.get(moi, moi)

            validity_entry = {
                "diseaseId": disease_id,
                "disease": disease_label,
                "inheritance": inheritance,
                "classification": classification,
                "classificationDate": classification_date,
            }

            if gene_symbol not in gene_to_validity:
                gene_to_validity[gene_symbol] = []
            gene_to_validity[gene_symbol].append(validity_entry)

    logging.info("Loaded ClinGen Gene-Disease Validity for %d genes", len(gene_to_validity))
    return gene_to_validity


def load_gnomad_gene_constraints(gnomad_constraints_path):
    """
    Load gnomAD gene constraint metrics (LOEUF, pLI, etc.) from TSV.

    For genes with multiple transcripts, keeps the entry with the lowest LOEUF
    (most constrained), matching the VEP LOEUF plugin's gene-level behavior.

    Returns a dict mapping gene symbol to constraint data matching Nirvana format:
        {pLi, pRec, pNull, synZ, misZ, loeuf}

    Args:
        gnomad_constraints_path (str): Path to gnomad.v2.1.1.lof_metrics.by_gene.txt

    Returns:
        dict: gene_symbol -> gnomAD constraint dict
    """
    gene_to_gnomad = {}

    def safe_float(val):
        try:
            return float(val)
        except (ValueError, TypeError):
            return None

    with open(gnomad_constraints_path, 'r') as f:
        header = f.readline().strip().split('\t')
        col = {name: i for i, name in enumerate(header)}

        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < len(header):
                continue

            gene = fields[col['gene']]
            loeuf = safe_float(fields[col['oe_lof_upper']])
            if loeuf is None:
                continue

            gnomad_entry = {
                "loeuf": loeuf,
                "pLi": safe_float(fields[col['pLI']]),
                "pRec": safe_float(fields[col['pRec']]),
                "pNull": safe_float(fields[col['pNull']]),
                "synZ": safe_float(fields[col['syn_z']]),
                "misZ": safe_float(fields[col['mis_z']]),
                "oe_mis_upper": safe_float(fields[col['oe_mis_upper']]),
            }

            # Keep the most constrained (lowest LOEUF) entry per gene
            existing = gene_to_gnomad.get(gene)
            if existing is None or loeuf < existing.get("loeuf", float('inf')):
                gene_to_gnomad[gene] = gnomad_entry

    logging.info("Loaded gnomAD constraint data for %d genes", len(gene_to_gnomad))
    return gene_to_gnomad


def build_hgnc_to_gene_data(clingen_csv_path, gnomad_constraints_path):
    """
    Build a complete hgnc_to_gene_data dict from external data sources.

    This is the annotator-independent alternative to loading gene data from
    Nirvana's embedded genes section.

    Args:
        clingen_csv_path (str): Path to ClinGen Gene-Disease Validity CSV
        gnomad_constraints_path (str): Path to gnomAD gene constraint TSV

    Returns:
        dict: gene_symbol -> {name, clingenGeneValidity, gnomAD}
    """
    gene_to_validity = load_clingen_gene_validity(clingen_csv_path)
    gene_to_gnomad = load_gnomad_gene_constraints(gnomad_constraints_path)

    # Merge into a single dict
    all_genes = set(gene_to_validity.keys()) | set(gene_to_gnomad.keys())
    hgnc_to_gene_data = {}

    for gene in all_genes:
        entry = {"name": gene}
        if gene in gene_to_validity:
            entry["clingenGeneValidity"] = gene_to_validity[gene]
        if gene in gene_to_gnomad:
            entry["gnomAD"] = gene_to_gnomad[gene]
        hgnc_to_gene_data[gene] = entry

    logging.info("Built gene data for %d genes (%d with ClinGen validity, %d with gnomAD constraints)",
                 len(hgnc_to_gene_data), len(gene_to_validity), len(gene_to_gnomad))
    return hgnc_to_gene_data
