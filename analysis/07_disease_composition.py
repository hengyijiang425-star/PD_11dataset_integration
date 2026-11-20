"""Composition of each Lewy-body diagnosis (and PD) vs control, per region.
Input  results/pseudobulk/pb_meta.csv (nuclei per donor x cell type; sorted libraries and Unresolved excluded)
Output results/followup/disease_vs_control_composition.csv
Model as 04_composition_test.py: empirical logit ~ disease + C(dataset), datasets holding both groups; BH within region x disease.
Run from the project root."""
import numpy as np, pandas as pd, statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

def region_group(r):
    r = str(r).lower()
    if any(k in r for k in ["nigra", "midbrain", "snpc"]): return "Midbrain"
    if any(k in r for k in ["cortex", "cingulate", "frontal"]): return "Cortex"
    return "Striatum/BF"

def run(pb_meta="results/pseudobulk/pb_meta.csv", min_donors=8):
    m = pd.read_csv(pb_meta); m["region_group"] = m.region.map(region_group)
    W = m.pivot_table(index=["region_group", "dataset", "donor", "diagnosis"], columns="cell_type", values="n_cells", aggfunc="sum", fill_value=0)
    n = W.sum(axis=1); P = W.div(n, axis=0).reset_index(); P["n"] = n.values
    rows = []
    for (reg, dx), _ in P[P.diagnosis != "Control"].groupby(["region_group", "diagnosis"]):
        s = P[(P.region_group == reg) & P.diagnosis.isin([dx, "Control"])].copy(); s["D"] = (s.diagnosis == dx).astype(int)
        s = s[s.dataset.isin([d for d, x in s.groupby("dataset") if x.D.nunique() == 2])]
        if s.D.sum() < 3: continue
        form = "y ~ D + C(dataset)" if s.dataset.nunique() > 1 else "y ~ D"
        for ct in W.columns:
            if (s[ct] > 0).sum() < min_donors: continue
            k = s[ct] * s.n; s["y"] = np.log((k + .5) / (s.n - k + .5)); f = smf.ols(form, s).fit()
            rows.append(dict(region=reg, disease=dx, cell_type=ct, n_dis=int(s.D.sum()), n_ctrl=int((1 - s.D).sum()),
                             datasets=", ".join(sorted(s.dataset.unique())), pct_dis=round(100 * s.loc[s.D == 1, ct].mean(), 2),
                             pct_ctrl=round(100 * s.loc[s.D == 0, ct].mean(), 2), coef=round(f.params["D"], 3), p=f.pvalues["D"]))
    R = pd.DataFrame(rows)
    R["padj"] = R.groupby(["region", "disease"]).p.transform(lambda p: multipletests(p, method="fdr_bh")[1])
    return R.sort_values(["region", "disease", "p"])

if __name__ == "__main__":
    import os; os.makedirs("results/followup", exist_ok=True)
    r = run(); r.to_csv("results/followup/disease_vs_control_composition.csv", index=False); print(r.head(20).to_string(index=False))
