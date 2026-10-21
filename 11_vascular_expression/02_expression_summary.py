"""Pseudobulk summary of a single-cell dataset for the locus genes.

For each donor x region x cell type: number of cells, mean log1p(CP10k), mean CP10k
and fraction of cells with any expression, for the genes in
results/vascular/genes_for_expression.tsv. Run from the repository root:

    python 11_vascular_expression/02_expression_summary.py hca_artery
    python 11_vascular_expression/02_expression_summary.py brain_aneurysm

Output: results/vascular/expr_{dataset}.tsv.gz
"""

import gzip
import os
import sys

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml


def load_config():
    with open("config.yml") as f:
        cfg = yaml.safe_load(f)
    if os.path.exists("config.local.yml"):
        with open("config.local.yml") as f:
            local = yaml.safe_load(f) or {}
        cfg["vascular"]["datasets"].update(local.get("vascular", {}).get("datasets", {}))
    return cfg


def read_elem(g):
    try:
        from anndata.io import read_elem as re_
    except ImportError:
        from anndata.experimental import read_elem as re_
    return re_(g)


def summarise(groups, n_genes, rows):
    """rows yields (cell_index_array, dense log1p block, dense cp10k block)."""
    codes, uniq = pd.factorize(groups)
    k = len(uniq)
    s_log = np.zeros((k, n_genes))
    s_lin = np.zeros((k, n_genes))
    s_nz = np.zeros((k, n_genes))
    for idx, lg, lin in rows:
        onehot = sp.csr_matrix((np.ones(len(idx)), (codes[idx], np.arange(len(idx)))), shape=(k, len(idx)))
        s_log += onehot @ lg
        s_lin += onehot @ lin
        s_nz += onehot @ (lin > 0)
    n = np.bincount(codes, minlength=k).astype(float)[:, None]
    return uniq, n[:, 0], s_log / n, s_lin / n, s_nz / n


def h5ad_rows(path, gene_ids, id_type, step=20000):
    f = h5py.File(path, "r")
    if "raw" in f and "X" in f["raw"]:
        X, var = f["raw/X"], read_elem(f["raw/var"])
    else:
        X, var = f["X"], read_elem(f["var"])
    obs = read_elem(f["obs"])
    names = var.index.astype(str) if id_type == "ensembl" else var["feature_name"].astype(str)
    names = pd.Index([n.split(".")[0] for n in names])
    col = names.get_indexer(gene_ids)
    found = col >= 0
    if not isinstance(X, h5py.Group) or X.attrs.get("encoding-type") != "csr_matrix":
        sys.exit("expected a CSR matrix in " + path)
    n_cells, n_var = X.attrs["shape"]
    indptr = X["indptr"][:]
    data, indices = X["data"], X["indices"]

    sample = data[: min(len(data), 100000)]
    if not np.allclose(sample, np.round(sample)):
        sys.exit("matrix in " + path + " does not look like raw counts")

    def rows():
        for s in range(0, n_cells, step):
            e = min(s + step, n_cells)
            a, b = indptr[s], indptr[e]
            m = sp.csr_matrix((data[a:b], indices[a:b], indptr[s:e + 1] - a), shape=(e - s, n_var))
            lib = np.asarray(m.sum(axis=1)).ravel()
            lib[lib == 0] = 1
            sub = m[:, col[found]].toarray() / lib[:, None] * 1e4
            lin = np.zeros((e - s, len(gene_ids)))
            lin[:, found] = sub
            yield np.arange(s, e), np.log1p(lin), lin

    return obs, rows(), found


def cellbrowser_rows(matrix, meta, genes, counts):
    obs = pd.read_csv(meta, sep="\t", index_col=0)
    with gzip.open(matrix, "rt") as f:
        cells = f.readline().rstrip("\n").split("\t")[1:]
        keep = {g: i for i, g in enumerate(genes)}
        block = np.zeros((len(cells), len(genes)), dtype=np.float32)
        lib = np.zeros(len(cells))
        found = np.zeros(len(genes), dtype=bool)
        for line in f:
            name, vals = line.split("\t", 1)
            sym = name.split("|")[-1]
            if counts or sym in keep:
                x = np.array(vals.split("\t"), dtype=float)
                if counts:
                    lib += x
                if sym in keep:
                    block[:, keep[sym]] = x
                    found[keep[sym]] = True
    obs = obs.reindex(cells)
    if counts:
        lib[lib == 0] = 1
        lin = block / lib[:, None] * 1e4
        lg = np.log1p(lin)
    else:
        # already log-normalised
        lg = block
        lin = np.expm1(block)
    return obs, iter([(np.arange(len(cells)), lg, lin)]), found


def mtx_rows(d, gene_ids):
    """10x-style matrix market (features x cells, raw counts) with a cell metadata table."""
    import scipy.io
    feats = pd.read_csv(d["features"], sep="\t", header=None)
    col_ix = 0 if d["gene_id"] == "ensembl" or feats.shape[1] == 1 else 1
    names = pd.Index([str(n).split(".")[0] if d["gene_id"] == "ensembl" else str(n)
                      for n in feats[col_ix]])
    cells = pd.read_csv(d["barcodes"], sep="\t", header=None)[0].astype(str).values
    with gzip.open(d["matrix"], "rb") as f:
        X = scipy.io.mmread(f).T.tocsr()          # cells x features
    if X.shape != (len(cells), len(names)):
        sys.exit(f"matrix is {X.shape}, expected {len(cells)} cells x {len(names)} features")
    if not np.allclose(X.data[:100000], np.round(X.data[:100000])):
        sys.exit("matrix in " + d["matrix"] + " does not look like raw counts")
    col = names.get_indexer(gene_ids)
    found = col >= 0
    lib = np.asarray(X.sum(axis=1)).ravel()
    lib[lib == 0] = 1
    lin = np.zeros((len(cells), len(gene_ids)))
    lin[:, found] = X[:, col[found]].toarray() / lib[:, None] * 1e4
    obs = pd.read_csv(d["meta"], sep="\t", index_col=0)
    missing = (~pd.Index(cells).isin(obs.index)).sum()
    if missing:
        print(f"{missing} of {len(cells)} cells have no metadata and are dropped")
    obs = obs.reindex(cells)
    keep = np.where(obs.notna().any(axis=1).values)[0]
    return obs.iloc[keep], iter([(np.arange(len(keep)), np.log1p(lin[keep]), lin[keep])]), found


def main():
    key = sys.argv[1]
    cfg = load_config()
    d = cfg["vascular"]["datasets"][key]
    out_dir = os.path.join(cfg["results_dir"], "vascular")
    genes = pd.read_csv(os.path.join(out_dir, "genes_for_expression.tsv"), sep="\t")

    if d["format"] == "h5ad":
        ids = genes["gene_id"] if d["gene_id"] == "ensembl" else genes["gene"]
        obs, rows, found = h5ad_rows(d["file"], list(ids), d["gene_id"])
    elif d["format"] == "mtx":
        ids = genes["gene_id"] if d["gene_id"] == "ensembl" else genes["gene"]
        obs, rows, found = mtx_rows(d, list(ids))
    else:
        obs, rows, found = cellbrowser_rows(d["matrix"], d["meta"], list(genes["gene"]), d.get("counts", False))

    for c in (d["donor"], d["region"], d["cell_type"]):
        if c not in obs.columns:
            sys.exit(f"column {c} not in the cell metadata; available: {', '.join(obs.columns)}")
    # several samples from one donor (e.g. ctrl099_21, ctrl099_22) are pooled by donor
    if d.get("donor_strip"):
        obs[d["donor"]] = obs[d["donor"]].astype(str).str.replace(d["donor_strip"], "", regex=True)
        print("donors:", ", ".join(sorted(obs[d["donor"]].dropna().unique())))
    print(f"{key}: {len(obs)} cells, {found.sum()} of {len(genes)} genes found")

    groups = (obs[d["donor"]].astype(str) + "\t" + obs[d["region"]].astype(str) + "\t" +
              obs[d["cell_type"]].astype(str)).values
    uniq, n, m_log, m_lin, pct = summarise(groups, len(genes), rows)

    lab = pd.DataFrame([u.split("\t") for u in uniq], columns=["donor", "region", "cell_type"])
    res = []
    for j in np.where(found)[0]:
        res.append(lab.assign(gene_id=genes["gene_id"][j], gene=genes["gene"][j], n_cells=n.astype(int),
                              mean_log=m_log[:, j], mean_cp10k=m_lin[:, j], pct=pct[:, j]))
    res = pd.concat(res)
    res.to_csv(os.path.join(out_dir, f"expr_{key}.tsv.gz"), sep="\t", index=False, float_format="%.5g")


if __name__ == "__main__":
    main()
