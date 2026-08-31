"""
Derives constraint-stratified AF cutoffs for ACMG evidence rules using
gene-normalized KDE, Tavtigian LR thresholds, and cluster bootstrap CIs.
Trains on a 90% gene-level split, evaluates on 10% holdout, and outputs
cutoff and holdout TSVs for downstream plotting (script 09).
"""

import argparse
import csv
import numpy as np

# Adaptive window minimum counts, scale down at edges of constraint range
MIN_PATH = 300
MIN_BENIGN = 300
MIN_GENES = 15

KDE_BW = 0.03
SMOOTH_BW = 0.10
BOOTSTRAP_ITER = 400
# Point estimates need more samples than bootstrap resamples (CIs tolerate noise)
MIN_WINDOW_SAMPLES = 10
MIN_BOOT_SAMPLES = 5
DEPLOY_BIN_WIDTH = 0.1
HOLDOUT_FRAC = 0.10
SEED = 42

N_BINS = 1000
AF_GRID = np.linspace(-7, 0, N_BINS)
BIN_WIDTH = AF_GRID[1] - AF_GRID[0]

# Tavtigian LR thresholds, benign rules use 1/threshold, pathogenic rules use threshold directly
benign_thresholds = {'BA1': 1124.0, 'BS1': 33.53, 'BS1_Moderate': 5.79, 'BS1_Supporting': 2.406}
path_thresholds = {'PM2_Supporting': 2.406}
all_rules = list(benign_thresholds.keys()) + list(path_thresholds.keys())

# Signed evidence points on the Tavtigian scale as log2 LR
evidence_points = {
    'BA1': -np.log2(1124.0),
    'BS1': -np.log2(33.53),
    'BS1_Moderate': -np.log2(5.79),
    'BS1_Supporting': -np.log2(2.406),
    'PM2_Supporting': +np.log2(2.406),
    'VUS': 0.0,
}

XLIM = (0.15, 1.8)


# Weighted KDE via histogram binning + Gaussian convolution
def evaluate_kde(data, weights, bandwidth):
    bin_indices = np.clip(((data - AF_GRID[0]) / BIN_WIDTH).astype(int), 0, N_BINS - 1)
    hist = np.zeros(N_BINS)
    np.add.at(hist, bin_indices, weights)
    kernel_half = int(4 * bandwidth / BIN_WIDTH)
    kernel_x = np.arange(-kernel_half, kernel_half + 1) * BIN_WIDTH
    kernel = np.exp(-0.5 * (kernel_x / bandwidth) ** 2) / (bandwidth * np.sqrt(2 * np.pi))
    return np.convolve(hist, kernel, mode='same') / BIN_WIDTH


# Per-variant weights so each gene contributes equally
def gene_weights(genes):
    counts = {}
    for g in genes:
        counts[g] = counts.get(g, 0) + 1
    n = len(counts)
    return np.array([1.0 / (counts[g] * n) for g in genes])


# Gaussian kernel smoother in log10(AF) space
def smooth_curve(x, y, bw):
    x = np.asarray(x, dtype=float)
    log_y = np.log10(np.asarray(y, dtype=float))
    smoothed = np.empty_like(log_y)
    for i in range(len(x)):
        w = np.exp(-0.5 * ((x - x[i]) / bw) ** 2)
        w /= w.sum()
        smoothed[i] = np.sum(w * log_y)
    return 10 ** smoothed


# Find the last grid index where LR exceeds threshold, scanning in direction
def find_lr_crossing(lr, meaningful, threshold, direction='high_to_low'):
    if direction == 'high_to_low':
        for j in range(len(lr) - 1, -1, -1):
            if meaningful[j] and lr[j] > threshold:
                return j
    else:
        last_above = -1
        for j in range(len(lr)):
            if meaningful[j] and lr[j] >= threshold:
                last_above = j
            elif meaningful[j] and lr[j] < threshold and last_above >= 0:
                break
        return last_above
    return -1


# Find AF where LR crosses each Tavtigian threshold
def compute_cutoffs(path_afs, benign_afs, path_w, benign_w):
    pdf_path = evaluate_kde(path_afs, path_w, KDE_BW)
    pdf_benign = evaluate_kde(benign_afs, benign_w, KDE_BW)

    density_floor = 1e-4 * max(np.max(pdf_path), np.max(pdf_benign))
    meaningful = (pdf_path + pdf_benign) > density_floor

    lr = np.full_like(AF_GRID, np.nan)
    lr[meaningful] = pdf_path[meaningful] / (pdf_benign[meaningful] + 1e-12)

    cutoffs = {}
    for rule, threshold in benign_thresholds.items():
        idx = find_lr_crossing(lr, meaningful, 1.0 / threshold, 'high_to_low')
        cutoffs[rule] = 10 ** AF_GRID[idx + 1] if 0 <= idx < len(AF_GRID) - 1 else None

    for rule, threshold in path_thresholds.items():
        idx = find_lr_crossing(lr, meaningful, threshold, 'low_to_high')
        cutoffs[rule] = 10 ** AF_GRID[idx] if idx >= 0 else None

    return cutoffs


# Extract and smooth a single rules cutoff series across sliding window centers
def extract_series(results, rule):
    centers, vals, los, his = [], [], [], []
    for r in results:
        v = r[5].get(rule)
        if v is not None:
            centers.append(r[0])
            vals.append(v)
            los.append(r[6].get(rule, v))
            his.append(r[7].get(rule, v))
    c, v, lo, hi = np.array(centers), np.array(vals), np.array(los), np.array(his)
    if len(c) > 3:
        v = smooth_curve(c, v, SMOOTH_BW)
        lo = smooth_curve(c, lo, SMOOTH_BW)
        hi = smooth_curve(c, hi, SMOOTH_BW)
    return c, v, lo, hi


# Gene-level cluster bootstrap for confidence intervals
def bootstrap_cutoffs(afs, is_path, genes):
    gene_to_idx = {}
    for i, g in enumerate(genes):
        gene_to_idx.setdefault(g, []).append(i)
    unique_genes = list(gene_to_idx.keys())

    boot_cutoffs = {rule: [] for rule in all_rules}
    for _ in range(BOOTSTRAP_ITER):
        boot_genes = np.random.choice(unique_genes, size=len(unique_genes), replace=True)
        boot_idx = np.concatenate([gene_to_idx[g] for g in boot_genes])
        b_afs, b_path, b_genes = afs[boot_idx], is_path[boot_idx], genes[boot_idx]
        if np.sum(b_path) < MIN_BOOT_SAMPLES or np.sum(~b_path) < MIN_BOOT_SAMPLES:
            continue
        bc = compute_cutoffs(b_afs[b_path], b_afs[~b_path],
                             gene_weights(b_genes[b_path]), gene_weights(b_genes[~b_path]))
        for rule in all_rules:
            if bc.get(rule) is not None:
                boot_cutoffs[rule].append(bc[rule])

    lo, hi = {}, {}
    for rule in all_rules:
        vals = boot_cutoffs[rule]
        if len(vals) >= 10:
            lo[rule] = np.percentile(vals, 2.5)
            hi[rule] = np.percentile(vals, 97.5)
        else:
            lo[rule] = hi[rule] = None

    return boot_cutoffs, lo, hi


def load_calibration_dataset(input_path):
    ad_rows = []
    ar_rows = []
    n_filtered_imputed = 0
    with open(input_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row['inheritance_mode'] not in ('AD', 'AR'):
                continue
            if row['clinical_significance'] == 'pathogenic' and float(row['gnomad_popmax_af']) > 0.10:
                continue
            if int(row['gnomad_all_ac']) == 0:
                n_filtered_imputed += 1
                continue
            if row['inheritance_mode'] == 'AD':
                ad_rows.append(row)
            else:
                ar_rows.append(row)

    print(f"Loaded {len(ad_rows)} AD, {len(ar_rows)} AR (filtered {n_filtered_imputed} imputed)")
    return ad_rows, ar_rows


# Slide adaptive window along constraint axis, compute cutoffs with bootstrap CIs
# Upper CI corrects for noisy point estimates in genes with few expected variants
def run_analysis(variants, moi_label):
    oe_vals = np.array([float(v['mis_oe_ci_upper']) for v in variants])
    is_path = np.array([v['clinical_significance'] == 'pathogenic' for v in variants])
    log_afs = np.array([np.log10(float(v['gnomad_popmax_af'])) for v in variants])
    genes = np.array([v['gene_symbol'] for v in variants])

    oe_max = float(max(oe_vals))
    centers = np.arange(XLIM[0], min(oe_max, XLIM[1]) + 0.005, 0.005)
    results = []

    for ci, center in enumerate(centers):
        radius = 0.02
        while True:
            low, high = center - radius, center + radius
            effective_width = min(high, oe_max) - max(low, 0.0)
            edge_scale = max(effective_width / (high - low), 0.1)

            in_window = (oe_vals >= low) & (oe_vals <= high)
            n_path_w = np.sum(in_window & is_path)
            n_benign_w = np.sum(in_window & ~is_path)
            window_genes = set(genes[in_window])
            if (n_path_w >= MIN_PATH * edge_scale and
                n_benign_w >= MIN_BENIGN * edge_scale and
                len(window_genes) >= MIN_GENES):
                break
            radius += 0.005
            if radius > 3.0:
                break

        w_genes = genes[in_window]
        w_afs = log_afs[in_window]
        w_path = is_path[in_window]
        n_path, n_benign = int(n_path_w), int(n_benign_w)
        n_genes = len(window_genes)

        if n_path < MIN_WINDOW_SAMPLES or n_benign < MIN_WINDOW_SAMPLES:
            continue

        cutoffs = compute_cutoffs(w_afs[w_path], w_afs[~w_path],
                                  gene_weights(w_genes[w_path]), gene_weights(w_genes[~w_path]))

        _, lo, hi = bootstrap_cutoffs(w_afs, w_path, w_genes)
        results.append((center, radius, n_path, n_benign, n_genes, cutoffs, lo, hi))

    print(f"  {moi_label} stratified: {len(results)} center points")
    return results


# Single global window, bootstrap median as point estimate
def train_flat_baseline(variants, moi_label):
    is_path = np.array([v['clinical_significance'] == 'pathogenic' for v in variants])
    log_afs = np.array([np.log10(float(v['gnomad_popmax_af'])) for v in variants])
    genes = np.array([v['gene_symbol'] for v in variants])

    cutoffs = compute_cutoffs(log_afs[is_path], log_afs[~is_path],
                              gene_weights(genes[is_path]), gene_weights(genes[~is_path]))

    boot_cutoffs, lo, hi = bootstrap_cutoffs(log_afs, is_path, genes)
    for rule in all_rules:
        vals = boot_cutoffs[rule]
        if len(vals) >= 10:
            cutoffs[rule] = np.percentile(vals, 50)

    return cutoffs, lo, hi


# Gene-level split ensuring no data leakage between train and test
def holdout_split(variants, holdout_frac, seed):
    rng = np.random.RandomState(seed)
    genes = sorted(set(v['gene_symbol'] for v in variants))
    rng.shuffle(genes)
    n_holdout = max(1, int(len(genes) * holdout_frac))
    holdout_genes = set(genes[:n_holdout])
    train = [v for v in variants if v['gene_symbol'] not in holdout_genes]
    test = [v for v in variants if v['gene_symbol'] in holdout_genes]
    return train, test, holdout_genes


# Downsample majority class so test has equal benign and pathogenic counts
def balance_test_set(test_variants, seed):
    rng = np.random.RandomState(seed)
    benign = [v for v in test_variants if v['clinical_significance'] == 'benign']
    path = [v for v in test_variants if v['clinical_significance'] == 'pathogenic']
    n = min(len(benign), len(path))
    if len(benign) > n:
        benign = list(rng.choice(benign, size=n, replace=False))
    if len(path) > n:
        path = list(rng.choice(path, size=n, replace=False))
    balanced = benign + path
    rng.shuffle(balanced)
    return balanced


# Unified classification, threshold_fn(rule) returns the AF threshold or None
def assign_code(af, threshold_fn):
    for rule in ['BA1', 'BS1', 'BS1_Moderate', 'BS1_Supporting']:
        t = threshold_fn(rule)
        if t is not None and af >= t:
            return rule
    t = threshold_fn('PM2_Supporting')
    if t is not None and af <= t:
        return 'PM2_Supporting'
    return 'VUS'


# Build threshold lookup functions for each model type
def make_stratified_lookup(interp, constraint_val):
    def lookup(rule):
        if rule not in interp:
            return None
        c, v = interp[rule]
        if constraint_val < c[0] or constraint_val > c[-1]:
            return None
        return float(np.interp(constraint_val, c, v))
    return lookup


def make_flat_lookup(cutoffs):
    return lambda rule: cutoffs.get(rule)


def make_binned_lookup_fn(binned_lookup, bin_edges, constraint_val):
    # Early return handles out-of-range, so bin_idx is guaranteed in bounds
    if constraint_val < bin_edges[0] or constraint_val >= bin_edges[-1]:
        return lambda rule: None
    bin_idx = int((constraint_val - bin_edges[0]) / (bin_edges[1] - bin_edges[0]))
    def lookup(rule):
        v = binned_lookup[rule][bin_idx]
        return None if np.isnan(v) else v
    return lookup


# Build interpolation arrays from stratified results
def build_cutoff_interpolators(results):
    interp = {}
    for rule in all_rules:
        c, v, _, _ = extract_series(results, rule)
        if len(c) > 0:
            interp[rule] = (c, v)
    return interp


# Discretize stratified curves into constraint bins, trimming sparse edges
def build_binned_lookup(results, variants, bin_width=None):
    if bin_width is None:
        bin_width = DEPLOY_BIN_WIDTH
    all_bin_edges = np.arange(XLIM[0], XLIM[1] + bin_width / 2, bin_width)
    n_all = len(all_bin_edges) - 1

    genes_per_bin = [set() for _ in range(n_all)]
    for v in variants:
        cv = float(v['mis_oe_ci_upper'])
        bi = int((cv - all_bin_edges[0]) / bin_width)
        bi = max(0, min(bi, n_all - 1))
        genes_per_bin[bi].add(v['gene_symbol'])
    gene_counts = [len(s) for s in genes_per_bin]

    first_valid = 0
    while first_valid < n_all and gene_counts[first_valid] < MIN_GENES:
        first_valid += 1
    last_valid = n_all - 1
    while last_valid >= 0 and gene_counts[last_valid] < MIN_GENES:
        last_valid -= 1
    if first_valid > last_valid:
        first_valid, last_valid = 0, n_all - 1

    bin_edges = all_bin_edges[first_valid:last_valid + 2]
    n_bins = len(bin_edges) - 1

    lookup = {}
    for rule in all_rules:
        c, v, _, _ = extract_series(results, rule)
        if len(c) == 0:
            lookup[rule] = np.full(n_bins, np.nan)
            continue
        binned = np.full(n_bins, np.nan)
        for i in range(n_bins):
            lo_edge, hi_edge = bin_edges[i], bin_edges[i + 1]
            mask = (c >= lo_edge) & (c < hi_edge)
            if mask.any():
                binned[i] = np.median(v[mask])
            elif c[0] <= (lo_edge + hi_edge) / 2 <= c[-1]:
                binned[i] = float(np.interp((lo_edge + hi_edge) / 2, c, v))
        lookup[rule] = binned

    return lookup, bin_edges


# Asymmetric hinge loss penalizing evidence in the wrong direction
def compute_loss(true_label, code):
    y = 1.0 if true_label == 'pathogenic' else -1.0
    s = evidence_points[code]
    return max(0.0, -y * s)


# Score test variants and print summary metrics
def evaluate_model(test_variants, classify_fn, model_name):
    correct = incorrect = vus_count = 0
    total_loss = 0.0
    losses = []

    for v in test_variants:
        true_label = v['clinical_significance']
        code = classify_fn(v)
        loss = compute_loss(true_label, code)
        total_loss += loss
        losses.append(loss)

        if code == 'VUS':
            vus_count += 1
        elif (true_label == 'benign' and code in ('BA1', 'BS1', 'BS1_Moderate', 'BS1_Supporting')) or \
             (true_label == 'pathogenic' and code == 'PM2_Supporting'):
            correct += 1
        else:
            incorrect += 1

    total_judged = correct + incorrect
    total_test = len(test_variants)
    accuracy = correct / total_judged if total_judged > 0 else 0
    mean_loss = total_loss / total_test if total_test > 0 else 0

    print(f"  {model_name:35s} acc={accuracy:.1%}  loss={mean_loss:.4f}  coverage={total_judged}/{total_test} ({total_judged/total_test:.1%})")

    return {'accuracy': accuracy, 'mean_loss': mean_loss, 'correct': correct,
            'incorrect': incorrect, 'vus_count': vus_count, 'losses': losses}


# Classify a variant using the appropriate per-MOI model
def make_classifier(model_type, per_moi_data, no_moi_cutoffs=None):
    def classify(v):
        af = float(v['gnomad_popmax_af'])
        cv = float(v['mis_oe_ci_upper'])
        moi = v['inheritance_mode']

        if model_type == 'stratified':
            return assign_code(af, make_stratified_lookup(per_moi_data[moi]['train_interp'], cv))
        elif model_type == 'binned':
            return assign_code(af, make_binned_lookup_fn(
                per_moi_data[moi]['binned_lookup'], per_moi_data[moi]['binned_edges'], cv))
        elif model_type == 'flat_permoi':
            return assign_code(af, make_flat_lookup(per_moi_data[moi]['flat_cutoffs']))
        elif model_type == 'flat_nomoi':
            return assign_code(af, make_flat_lookup(no_moi_cutoffs))
    return classify


def write_cutoffs_tsv(path, per_moi_data, no_moi_cutoffs):
    fieldnames = ['model_type', 'moi', 'rule', 'center', 'cutoff_value',
                  'ci_lo', 'ci_hi', 'bin_lo', 'bin_hi']
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, delimiter='\t', fieldnames=fieldnames)
        w.writeheader()

        for moi_label in ('AD', 'AR'):
            data = per_moi_data[moi_label]

            for rule in all_rules:
                c, v, lo_arr, hi_arr = extract_series(data['train_results'], rule)
                for i in range(len(c)):
                    w.writerow({
                        'model_type': 'stratified', 'moi': moi_label, 'rule': rule,
                        'center': f'{c[i]:.4f}', 'cutoff_value': f'{v[i]:.6e}',
                        'ci_lo': f'{lo_arr[i]:.6e}', 'ci_hi': f'{hi_arr[i]:.6e}',
                        'bin_lo': '', 'bin_hi': '',
                    })

            bl = data['binned_lookup']
            be = data['binned_edges']
            for rule in all_rules:
                vals = bl[rule]
                for bi in range(len(be) - 1):
                    if not np.isnan(vals[bi]):
                        w.writerow({
                            'model_type': 'binned', 'moi': moi_label, 'rule': rule,
                            'center': '', 'cutoff_value': f'{vals[bi]:.6e}',
                            'ci_lo': '', 'ci_hi': '',
                            'bin_lo': f'{be[bi]:.4f}', 'bin_hi': f'{be[bi+1]:.4f}',
                        })

            fc = data['flat_cutoffs']
            for rule in all_rules:
                v = fc.get(rule)
                if v is not None:
                    w.writerow({
                        'model_type': 'flat', 'moi': moi_label, 'rule': rule,
                        'center': '', 'cutoff_value': f'{v:.6e}',
                        'ci_lo': '', 'ci_hi': '', 'bin_lo': '', 'bin_hi': '',
                    })

        for rule in all_rules:
            v = no_moi_cutoffs.get(rule)
            if v is not None:
                w.writerow({
                    'model_type': 'flat', 'moi': 'ALL', 'rule': rule,
                    'center': '', 'cutoff_value': f'{v:.6e}',
                    'ci_lo': '', 'ci_hi': '', 'bin_lo': '', 'bin_hi': '',
                })

    print(f"Wrote {path}")


def write_holdout_tsv(path, all_test, classifiers):
    fieldnames = [
        'chromosome', 'position', 'ref_allele', 'alt_allele',
        'gene_symbol', 'inheritance_mode', 'clinical_significance',
        'gnomad_popmax_af', 'constraint_value',
        'code_stratified', 'code_binned', 'code_flat_permoi', 'code_flat_nomoi',
    ]
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, delimiter='\t', fieldnames=fieldnames)
        w.writeheader()
        for v in all_test:
            w.writerow({
                'chromosome': v['chromosome'],
                'position': v['position'],
                'ref_allele': v['ref_allele'],
                'alt_allele': v['alt_allele'],
                'gene_symbol': v['gene_symbol'],
                'inheritance_mode': v['inheritance_mode'],
                'clinical_significance': v['clinical_significance'],
                'gnomad_popmax_af': v['gnomad_popmax_af'],
                'constraint_value': v['mis_oe_ci_upper'],
                'code_stratified': classifiers['stratified'](v),
                'code_binned': classifiers['binned'](v),
                'code_flat_permoi': classifiers['flat_permoi'](v),
                'code_flat_nomoi': classifiers['flat_nomoi'](v),
            })

    print(f"Wrote {path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Path to calibration_dataset TSV")
    parser.add_argument("--output_dir", required=True, help="Directory for output plots and TSVs")
    args = parser.parse_args()

    np.random.seed(SEED)
    ad_rows, ar_rows = load_calibration_dataset(args.input)

    per_moi_data = {}

    for moi_label, rows in [('AD', ad_rows), ('AR', ar_rows)]:
        train, test_raw, holdout_genes = holdout_split(rows, HOLDOUT_FRAC, SEED)
        test = balance_test_set(test_raw, SEED)

        train_results = run_analysis(train, moi_label)
        train_interp = build_cutoff_interpolators(train_results)
        flat_cutoffs, _, _ = train_flat_baseline(train, f"{moi_label}-flat")
        binned_lookup, binned_edges = build_binned_lookup(train_results, train)

        per_moi_data[moi_label] = {
            'train_results': train_results, 'train_interp': train_interp,
            'flat_cutoffs': flat_cutoffs,
            'binned_lookup': binned_lookup, 'binned_edges': binned_edges,
            'train': train, 'test': test,
        }

        for mt, label in [('stratified', 'Stratified'), ('binned', f'Binned ({DEPLOY_BIN_WIDTH}-wide)'), ('flat_permoi', 'Flat baseline')]:
            evaluate_model(test, make_classifier(mt, per_moi_data), f"{moi_label} {label}")

    all_train = []
    for moi_label in ('AD', 'AR'):
        all_train.extend(per_moi_data[moi_label]['train'])
    no_moi_cutoffs, _, _ = train_flat_baseline(all_train, "no-MOI-flat")

    all_test = []
    for moi_label in ('AD', 'AR'):
        all_test.extend(per_moi_data[moi_label]['test'])
    classifiers = {mt: make_classifier(mt, per_moi_data, no_moi_cutoffs)
                   for mt in ('stratified', 'binned', 'flat_permoi', 'flat_nomoi')}

    model_labels = {
        'stratified': 'Per-MOI + stratified',
        'binned': f'Per-MOI + binned ({DEPLOY_BIN_WIDTH}-wide)',
        'flat_permoi': 'Per-MOI flat baseline',
        'flat_nomoi': 'No-MOI flat baseline',
    }
    for mt, label in model_labels.items():
        evaluate_model(all_test, classifiers[mt], label)

    write_cutoffs_tsv(f"{args.output_dir}/mis_oe_upper_cutoffs.tsv", per_moi_data, no_moi_cutoffs)
    write_holdout_tsv(f"{args.output_dir}/mis_oe_upper_holdout_variants.tsv", all_test, classifiers)

    print("Done")


if __name__ == '__main__':
    main()
