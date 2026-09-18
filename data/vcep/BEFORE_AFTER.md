# VCEP integration — before / after F1

Both runs use the same inputs (Nirvana annotation of the eRepo dataset,
10,164 variants, `test/data/erepo_03.02.2026.vcf`). The only difference is
whether `BA1_BS1_PM2_vcep_af_cutoffs_fp` is present in the paths JSON.

Regenerate with:

    python3 bias_2015.py test/data/erepo_02.09.2026.nirvana.json <paths.json> <out.tsv>
    python3 src/scripts/evaluate_bias_erepo_validation.py test/data/erepo_02.09.2026.tsv <out.tsv> test/data/erepo_03.02.2026.vcf bias --sens f1 --conc f2 --full_cod f3

## The three targeted codes

| Code | Metric    | Baseline (LOEUF-only) | With VCEP | Δ         |
|------|-----------|-----------------------|-----------|-----------|
| BA1  | TP        | 511                   | 718       | +207      |
| BA1  | FP        | 50                    | 183       | +133      |
| BA1  | FN        | 442                   | 235       | −207      |
| BA1  | Precision | 0.9109                | 0.7969    | −0.114    |
| BA1  | Recall    | 0.5362                | 0.7534    | **+0.217** |
| BA1  | **F1**    | 0.6750                | **0.7745** | **+0.100** |
| BS1  | TP        | 292                   | 253       | −39       |
| BS1  | FP        | 571                   | 359       | −212      |
| BS1  | FN        | 353                   | 392       | +39       |
| BS1  | Precision | 0.3384                | 0.4134    | +0.075    |
| BS1  | Recall    | 0.4527                | 0.3922    | −0.061    |
| BS1  | **F1**    | 0.3873                | **0.4025** | **+0.015** |
| PM2  | TP        | 6000                  | 6396      | +396      |
| PM2  | FP        | 1033                  | 719       | −314      |
| PM2  | FN        | 782                   | 386       | −396      |
| PM2  | Precision | 0.8531                | 0.8989    | +0.046    |
| PM2  | Recall    | 0.8847                | 0.9431    | +0.058    |
| PM2  | **F1**    | 0.8686                | **0.9205** | **+0.052** |

## Whole-pipeline effects

| Metric                        | Baseline | With VCEP | Δ         |
|-------------------------------|----------|-----------|-----------|
| Total F1 across 28 codes      | 8.9007   | **9.0673** | **+0.167** |
| Correctly-called codes        | 17,022   | 17,586    | +564      |
| Codes missed (FN)             | 14,098   | 13,534    | −564      |
| Codes called incorrectly (FP) | 17,567   | 17,174    | −393      |
| Pathogenic sensitivity        | 0.8943   | 0.8855    | −0.009    |
| Pathogenic specificity        | 0.8795   | 0.8836    | +0.004    |
| Benign sensitivity            | 0.8768   | 0.9029    | +0.026    |
| Benign specificity            | 0.9712   | 0.9632    | −0.008    |

## Codes not touched by this change (sanity)

Unchanged codes — should be byte-identical between the two runs; confirms VCEP
isn't leaking into unrelated classifiers:

| Code | Before F1 | After F1 |
|------|-----------|----------|
| PVS1 | 0.8975    | 0.8975   |
| PS3  | 0.4206    | 0.4206   |

## Interpretation

- **PM2 is the biggest win.** F1 +0.052 (0.87 → 0.92) with both precision and
  recall improving. Comes primarily from the "absent from controls" semantic
  rule (49 previously-unparsed VCEP entries now fire on gnomAD-absent variants)
  and from popmax replacing raw AF on genes with FAF-based rules.
- **BA1 recall jumps ~+22 points** (0.54 → 0.75) at the cost of ~11 points of
  precision. Net F1 improves by +0.10. VCEPs use lower BA1 thresholds than the
  LOEUF-tiered defaults on many disease genes, so we call more benign variants
  correctly. The precision cost is real but small in absolute terms (+133 FPs
  against +207 recovered TPs).
- **BS1 is the smallest and noisiest change** (+0.015 F1). Most VCEPs specify
  BS1 in terms of FAF95, which BIAS approximates with popmax. That's the right
  metric family but a slight over-estimate, hence the modest gain.
- **Benign sensitivity +2.6 points, pathogenic sensitivity −0.9 points** — the
  direction matches the design intent (VCEPs let us call more variants benign
  in covered genes). Specificity moves are negligible in both directions.
- **PVS1 and PS3 are byte-identical**, confirming the change is scoped to the
  three AF-based codes as intended.

## Coverage of the input table

- CSpec sequence-variant-interpretation docs surveyed: **~205** (index)
- Released / approved docs actually parsed: **~118**
- Unique genes with at least one BA1/BS1/PM2 rule: **130**
- Rules emitted to the TSV: **451**
- Rules in the review sidecar (unparseable): **8**
- **Parse success rate: 98.3 %**
