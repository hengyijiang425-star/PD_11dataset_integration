"""PD vs control cell-type composition test, per region.

Input  results/composition_test/donor_celltype_proportions.csv  (one row per donor; sorted /
       neuron-enriched libraries already removed), results/pseudobulk/pb_meta.csv (nuclei per donor)
Output results/composition_test/PD_vs_control_composition_test.csv

Model, per region and cell type:  empirical logit of the donor's proportion,
    y = log((k + 0.5) / (n - k + 0.5)),  k = nuclei of that type, n = annotated nuclei of the donor
    y ~ PD + C(dataset)   (OLS; dataset dropped when the region has one dataset)
BH-FDR within region. Only cell types present (>0) in at least 10 donors of the region are tested.
Run from the project root: python analysis/04_composition_test.py
"""
import numpy as np, pandas as pd, statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

def run(props="results/composition_test/donor_celltype_proportions.csv",
        pb_meta="results/pseudobulk/pb_meta.csv", min_donors=10):
    d = pd.read_csv(props)
    d = d[d.diagnosis.isin(["PD", "Control"])].copy()
    n = pd.read_csv(pb_meta).groupby("donor").n_cells.sum()
    d["n"] = d.donor.map(n)
    cts = [c for c in d.columns if c not in ("region_group", "dataset", "donor", "diagnosis", "n")]
    rows = []
    for reg, g in d.groupby("region_group"):
        g = g.dropna(subset=["n"]).copy()
        g["PD"] = (g.diagnosis == "PD").astype(int)
        if g.PD.nunique() < 2:
            continue
        form = "y ~ PD + C(dataset)" if g.dataset.nunique() > 1 else "y ~ PD"
        for ct in cts:
            if (g[ct] > 0).sum() < min_donors:
                continue
            k = g[ct] * g.n
            g["y"] = np.log((k + 0.5) / (g.n - k + 0.5))
            fit = smf.ols(form, g).fit()
            per_ds = g.groupby(["dataset", "PD"])[ct].mean().unstack()
            per_ds = per_ds.dropna()
            rows.append(dict(region=reg, cell_type=ct,
                             donors_PD=int(g.PD.sum()), donors_Ctrl=int((1 - g.PD).sum()),
                             pct_PD=round(100 * g.loc[g.PD == 1, ct].mean(), 2),
                             pct_Ctrl=round(100 * g.loc[g.PD == 0, ct].mean(), 2),
                             logit_diff_PD_vs_Ctrl=round(fit.params["PD"], 3), p=fit.pvalues["PD"],
                             datasets_lower_in_PD=int((per_ds[1] < per_ds[0]).sum()),
                             datasets_higher_in_PD=int((per_ds[1] > per_ds[0]).sum())))
    res = pd.DataFrame(rows)
    res["padj"] = res.groupby("region").p.transform(lambda p: multipletests(p, method="fdr_bh")[1])
    return res.sort_values(["region", "p"])

if __name__ == "__main__":
    r = run()
    r.to_csv("results/composition_test/PD_vs_control_composition_test.csv", index=False)
    print(r.to_string(index=False))
