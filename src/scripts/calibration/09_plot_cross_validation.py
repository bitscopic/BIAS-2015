"""
Plot 4-model 2x2 scatter comparison from pre-computed cutoffs and holdout
variant TSVs produced by script 08. All codes are read directly from the
holdout TSV, no model training or classification happens here.
"""

import argparse
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("--cutoffs_tsv", required=True, help="Cutoffs TSV from script 08")
parser.add_argument("--holdout_tsv", required=True, help="Holdout variants TSV from script 08")
parser.add_argument("--output_dir", required=True, help="Directory for output plots")
args = parser.parse_args()

AXIS_LABEL = "Missense O/E upper CI"
XLIM = (0.15, 1.8)
PREFIX = "mis_oe_upper"

rule_styles = {
    'BA1': ('#E91E63', '-', 2.5),
    'BS1': ('#9C27B0', '-', 2.2),
    'BS1_Moderate': ('#1565C0', '-', 2.2),
    'BS1_Supporting': ('#42A5F5', '-', 1.8),
    'PM2_Supporting': ('#E65100', '--', 2.2),
}

evidence_points = {
    'BA1': -np.log2(1124.0),
    'BS1': -np.log2(33.53),
    'BS1_Moderate': -np.log2(5.79),
    'BS1_Supporting': -np.log2(2.406),
    'PM2_Supporting': +np.log2(2.406),
    'VUS': 0.0,
}


def load_cutoffs(path):
    # Load cutoffs TSV into nested dicts keyed by (model_type, moi, rule)
    stratified = defaultdict(lambda: defaultdict(lambda: {'centers': [], 'values': [], 'ci_lo': [], 'ci_hi': []}))
    binned = defaultdict(lambda: defaultdict(lambda: {'bin_lo': [], 'bin_hi': [], 'values': []}))
    flat = defaultdict(dict)

    with open(path, 'r') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            model_type = row['model_type']
            moi = row['moi']
            rule = row['rule']
            val = float(row['cutoff_value'])

            if model_type == 'stratified':
                s = stratified[moi][rule]
                s['centers'].append(float(row['center']))
                s['values'].append(val)
                s['ci_lo'].append(float(row['ci_lo']))
                s['ci_hi'].append(float(row['ci_hi']))
            elif model_type == 'binned':
                b = binned[moi][rule]
                b['bin_lo'].append(float(row['bin_lo']))
                b['bin_hi'].append(float(row['bin_hi']))
                b['values'].append(val)
            elif model_type == 'flat':
                flat[moi][rule] = val

    # Convert to arrays
    for moi in stratified:
        for rule in stratified[moi]:
            for key in ('centers', 'values', 'ci_lo', 'ci_hi'):
                stratified[moi][rule][key] = np.array(stratified[moi][rule][key])
    for moi in binned:
        for rule in binned[moi]:
            for key in ('bin_lo', 'bin_hi', 'values'):
                binned[moi][rule][key] = np.array(binned[moi][rule][key])

    return dict(stratified), dict(binned), dict(flat)


def load_holdout(path):
    # Load holdout variants TSV
    with open(path, 'r') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def compute_loss(true_label, code):
    y = 1.0 if true_label == 'pathogenic' else -1.0
    s = evidence_points[code]
    return max(0.0, -y * s)


def compute_panel_stats(variants, code_column):
    # Accuracy, loss, and coverage for one model
    correct = incorrect = vus_count = 0
    total_loss = 0.0

    for v in variants:
        code = v[code_column]
        true_label = v['clinical_significance']
        loss = compute_loss(true_label, code)
        total_loss += loss

        if code == 'VUS':
            vus_count += 1
        elif (true_label == 'benign' and code in ('BA1', 'BS1', 'BS1_Moderate', 'BS1_Supporting')) or \
             (true_label == 'pathogenic' and code == 'PM2_Supporting'):
            correct += 1
        else:
            incorrect += 1

    judged = correct + incorrect
    total = len(variants)
    return {
        'accuracy': correct / judged if judged > 0 else 0,
        'mean_loss': total_loss / total if total > 0 else 0,
        'judged': judged,
        'total': total,
    }


def plot_4model_scatter(variants, stratified, binned, flat, output_dir):
    # 2x2 scatter comparing all four models
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Helvetica Neue', 'Arial', 'DejaVu Sans'],
        'axes.spines.top': False,
        'axes.spines.right': False,
    })

    fig, ((ax_s, ax_bn), (ax_mf, ax_nf)) = plt.subplots(2, 2, figsize=(14, 14),
                                                          sharey=True, sharex=True)
    fig.patch.set_facecolor('white')

    code_columns = ['code_stratified', 'code_binned', 'code_flat_permoi', 'code_flat_nomoi']
    stats = {col: compute_panel_stats(variants, col) for col in code_columns}

    panel_labels = ['(a)', '(b)', '(c)', '(d)']
    panel_configs = [
        (ax_s, 'Per-MOI + constraint-stratified', 'code_stratified'),
        (ax_bn, 'Per-MOI + binned', 'code_binned'),
        (ax_mf, 'Per-MOI flat', 'code_flat_permoi'),
        (ax_nf, 'No-MOI flat (AD+AR pooled)', 'code_flat_nomoi'),
    ]

    for pi, (ax, model_name, code_col) in enumerate(panel_configs):
        ax.set_facecolor('white')

        # Draw cutoff lines or curves
        if code_col == 'code_stratified':
            for moi, ls_mod in [('AD', '-'), ('AR', '--')]:
                if moi not in stratified:
                    continue
                for rule, (color, _, lw) in rule_styles.items():
                    if rule not in stratified[moi]:
                        continue
                    s = stratified[moi][rule]
                    ax.plot(s['centers'], s['values'], color=color, linestyle=ls_mod,
                            linewidth=lw * 0.8, alpha=0.8 if moi == 'AD' else 0.5,
                            zorder=4, path_effects=[pe.withStroke(linewidth=lw + 1, foreground='white')])

        elif code_col == 'code_binned':
            for moi, ls_mod in [('AD', '-'), ('AR', '--')]:
                if moi not in binned:
                    continue
                for rule, (color, _, lw) in rule_styles.items():
                    if rule not in binned[moi]:
                        continue
                    b = binned[moi][rule]
                    for bi in range(len(b['values'])):
                        ax.plot([b['bin_lo'][bi], b['bin_hi'][bi]],
                                [b['values'][bi], b['values'][bi]],
                                color=color, linestyle=ls_mod, linewidth=lw * 0.8,
                                alpha=0.8 if moi == 'AD' else 0.5, zorder=4,
                                path_effects=[pe.withStroke(linewidth=lw + 1, foreground='white')])

        elif code_col == 'code_flat_permoi':
            for moi, ls_mod in [('AD', '-'), ('AR', '--')]:
                if moi not in flat:
                    continue
                for rule, (color, _, lw) in rule_styles.items():
                    v = flat[moi].get(rule)
                    if v is not None:
                        ax.axhline(v, color=color, linestyle=ls_mod, linewidth=lw * 0.8,
                                   alpha=0.8 if moi == 'AD' else 0.5, zorder=4,
                                   path_effects=[pe.withStroke(linewidth=lw + 1, foreground='white')])

        else:  # no-MOI flat
            if 'ALL' in flat:
                for rule, (color, ls, lw) in rule_styles.items():
                    v = flat['ALL'].get(rule)
                    if v is not None:
                        ax.axhline(v, color=color, linestyle=ls, linewidth=lw,
                                   zorder=4, path_effects=[pe.withStroke(linewidth=lw + 1.5, foreground='white')])

        # Scatter holdout variants
        for v in variants:
            af = float(v['gnomad_popmax_af'])
            constraint_val = float(v['constraint_value'])
            true_label = v['clinical_significance']
            code = v[code_col]

            if code == 'VUS':
                c, m, s, a = '#BDBDBD', 'o', 8, 0.12
            elif (true_label == 'benign' and code in ('BA1', 'BS1', 'BS1_Moderate', 'BS1_Supporting')) or \
                 (true_label == 'pathogenic' and code == 'PM2_Supporting'):
                c, m, s, a = '#4CAF50', 'o', 14, 0.25
            else:
                c, m, s, a = '#F44336', 'X', 30, 0.75
            ax.scatter(constraint_val, af, c=c, marker=m, s=s, alpha=a, zorder=3,
                       linewidths=0.3, edgecolors='white')

        ax.set_yscale('log')
        ax.set_xlim(*XLIM)
        ax.tick_params(axis='both', labelsize=11)
        ax.grid(True, which='major', alpha=0.15, linewidth=0.6, color='#666666')
        ax.grid(True, which='minor', alpha=0.06, linewidth=0.4, color='#999999')
        for spine in ax.spines.values():
            spine.set_color('#333333')
            spine.set_linewidth(0.6)

        st = stats[code_col]
        ax.set_title(
            f"{panel_labels[pi]}  {model_name}\n"
            f"Acc {st['accuracy']:.1%}  |  Loss {st['mean_loss']:.4f} bits/var  |  "
            f"Coverage {st['judged']}/{st['total']} ({st['judged']/st['total']:.1%})",
            fontsize=10, color='#333333', pad=8)

    ax_mf.set_xlabel(AXIS_LABEL, fontsize=13, labelpad=6)
    ax_nf.set_xlabel(AXIS_LABEL, fontsize=13, labelpad=6)
    ax_s.set_ylabel('gnomAD popmax AF', fontsize=13, labelpad=6)
    ax_mf.set_ylabel('gnomAD popmax AF', fontsize=13, labelpad=6)

    legend_lines = [
        Line2D([0], [0], color='#4CAF50', marker='o', linestyle='None', markersize=6, label='Correct'),
        Line2D([0], [0], color='#F44336', marker='X', linestyle='None', markersize=7, label='Incorrect'),
        Line2D([0], [0], color='#BDBDBD', marker='o', linestyle='None', markersize=5, label='VUS'),
        Line2D([0], [0], color='grey', linestyle='-', linewidth=1.5, label='AD'),
        Line2D([0], [0], color='grey', linestyle='--', linewidth=1.5, label='AR'),
    ]
    ax_bn.legend(handles=legend_lines, loc='lower left', fontsize=9.5, framealpha=0.95,
                 edgecolor='#CCCCCC', fancybox=False, borderpad=0.6, handletextpad=0.4)

    fig.subplots_adjust(hspace=0.22, wspace=0.08)
    scatter_path = f"{output_dir}/{PREFIX}_4model_scatter_combined.png"
    plt.savefig(scatter_path, dpi=300, facecolor='white', bbox_inches='tight')
    plt.close()
    print(f"Wrote {scatter_path}")


def main():
    stratified, binned, flat = load_cutoffs(args.cutoffs_tsv)
    variants = load_holdout(args.holdout_tsv)
    print(f"Loaded {len(variants)} holdout variants")

    n_mois = len(set(v['inheritance_mode'] for v in variants))
    n_benign = sum(1 for v in variants if v['clinical_significance'] == 'benign')
    n_path = sum(1 for v in variants if v['clinical_significance'] == 'pathogenic')
    print(f"  {n_benign} benign, {n_path} pathogenic, {n_mois} MOI groups")

    # Per-model summary
    for label, col in [('Stratified', 'code_stratified'), ('Binned', 'code_binned'),
                       ('Per-MOI flat', 'code_flat_permoi'), ('No-MOI flat', 'code_flat_nomoi')]:
        st = compute_panel_stats(variants, col)
        print(f"  {label:20s}  acc={st['accuracy']:.1%}  loss={st['mean_loss']:.4f}  "
              f"coverage={st['judged']}/{st['total']} ({st['judged']/st['total']:.1%})")

    plot_4model_scatter(variants, stratified, binned, flat, args.output_dir)
    print("Done.")


if __name__ == '__main__':
    main()
