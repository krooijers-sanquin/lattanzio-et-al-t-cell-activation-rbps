#!/usr/bin/env python
# coding: utf-8

from collections import defaultdict

import numpy as np
import pandas as pd
import duckdb
import blitzgsea
import xarray

np.random.seed(42)


def sluggify(s):
    from string import ascii_letters, digits
    return "".join(
        c for c in s if c in ascii_letters or c in digits or c in "-_."
    ).replace(":", "_").replace("-", "_").replace(".", "_")


sample_subset = pd.read_csv("pipeline_activation/samplesheet.tsv", sep="\t")

intensities = pd.read_csv(
    "pipeline_activation/intensity-values.tsv",
    sep="\t", index_col=0,
)

w_valid = intensities['w_valid']

hgnc_metadata = pd.read_csv(
    "pipeline_activation/raw/hgnc_metadata.tsv.gz",
    sep="\t", index_col=0,
)

# Load Uniprot -> HGNC mapping
DBFN = "data/vendor/KR20231013.human_omni_gene.db"
dbconn = duckdb.connect(DBFN, read_only=True)
uniprot2hgncid = (
    dbconn.execute("SELECT primary_accession, hgnc_refs[1] FROM uniprot WHERE length(hgnc_refs) == 1")
    .df()
    .set_index("primary_accession", verify_integrity=True)
    .iloc[:, 0]
    .to_dict()
)

# ### Load experimental data:
map_ds = xarray.load_dataset("pipeline_activation/model-output/global_map.h5")
E = xarray.concat(
    [map_ds['E'], map_ds['Ederiv'].rename({"deriv_coef": "coef"})],
    dim='coef',
)
w_expr = (
    (E.sel(coef='Intercept') >= -0.75)
    & (E.sel(coef='OOPS_abundance') >= -0.25)
)
assert w_expr.mean().item() >= 0.6

gene_symbols = hgnc_metadata.loc[E.coords['gene'].values]['symbols'].values
assert len(gene_symbols) == w_valid.sum()

coefslugs = defaultdict(lambda: set())
for coef in map(str, E.coords['coef'].values):
    slug = sluggify(coef)
    coefslugs[slug].add(coef)

assert all(len(v) == 1 for v in coefslugs.values())
coefslugs = {k: next(iter(v)) for (k, v) in coefslugs.items()}

available_collections = dbconn.execute(
    """
        SELECT gs.collection_name, COUNT(*) as n
        FROM msigdb_gene_set as gs
        GROUP BY gs.collection_name
        ORDER BY n
    """
    ).df()['collection_name'].values
# sluggify and make mapping:
collection_slugs = defaultdict(lambda: set())
for s in available_collections:
    slug = s.replace(":", "_")
    collection_slugs[slug].add(s)
assert all(len(v) == 1 for v in collection_slugs.values())
collection_slugs = {k: next(iter(v)) for (k, v) in collection_slugs.items()}

if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--outfile", required=True)
    ap.add_argument("-s", "--collection", required=True)
    ap.add_argument("-c", "--coef", required=True)

    args = ap.parse_args()

    assert args.coef in coefslugs
    real_coef = coefslugs[args.coef]

    assert args.collection in collection_slugs
    real_collection = collection_slugs[args.collection]

    gs_library = {
        k: list(set(v))
        for (k, v) in dbconn.execute(
            """
                SELECT gs.standard_name, s.symbol
                FROM msigdb_gene_set as gs
                INNER JOIN msigdb_gene_set_gene_symbol as gsgs ON gsgs.gene_set_id = gs.id
                INNER JOIN msigdb_gene_symbol as s ON gsgs.gene_symbol_id = s.id
                WHERE gs.collection_name = $collection
            """,
            {"collection": real_collection},
        ).df().groupby("standard_name")['symbol']
        if len(set(v)) >= 5
    }

    signature_df = pd.DataFrame(
        {
            0: gene_symbols,
            1: E.sel(coef=real_coef),
        },
    )[w_expr.values]

    # Run blitzgsea:
    try:
        result = blitzgsea.gsea(
            signature_df, gs_library,
            processes=4, progress=True, verbose=True,
            signature_cache=False,
            permutations=10_000, seed=42,
            add_noise=True,
            symmetric=False, anchors=50,
            kl_bins=200,
        )

        (
            result
            .assign(**{"abs_nes": lambda df: np.abs(df['nes'])})
            .sort_values("abs_nes", ascending=False)
        ).to_csv(args.outfile, sep="\t", index=True, header=True)
    except Exception:
        with open(args.outfile, 'wt') as fh:
            fh.write("Error\n")
