#!/usr/bin/env python
# coding: utf-8


import re
from collections import defaultdict
from functools import reduce
import logging
from pathlib import Path

import pandas as pd
import networkx as nx
import duckdb

# from tqdm.auto import tqdm
from rich.progress import track
from rich.logging import RichHandler
import warnings
from tqdm import TqdmExperimentalWarning


warnings.filterwarnings("ignore", category=TqdmExperimentalWarning)


# What the script _should_ do
# - Take a set of files with peptide level data, each possibly describing *multiple* samples
# - Align the peptides to possible proteins and assign isoforms as unambiguous as possible
# - Write wide-format table
# - (PERHAPS MOVE TO OTHER SCRIPT) Take the supplied annotations (by proteomics facility software)

def tqdm(*args, **kwargs):
    desc = kwargs.pop("desc", None)
    _ = kwargs.pop("unit", None)
    return track(*args, **kwargs, transient=True, description=desc)


def window_overlapping(s, w, n):
    slen = len(s)
    for i in range(0, slen - w + 1, n):
        yield s[i:(i + w)]


def natural_sort_key(s, _nsre=re.compile('([0-9]+)')):
    return [int(text) if text.isdigit() else text.lower()
            for text in _nsre.split(s)]


# Set up logging
FORMAT = "%(message)s"
logging.basicConfig(
    level="NOTSET", format=FORMAT, datefmt="[%X]", handlers=[RichHandler()]
)

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--outbase", required=True)
    ap.add_argument("-d", "--db", metavar="OMNI_GENE_DUCKDB", required=True)
    ap.add_argument("infile", nargs="+")

    args = ap.parse_args()

    outbase = Path(args.outbase)

    fns = args.infile
    assert fns
    log.info("Got files: %s" % ", ".join(fns))

    # Read Uniprot metadata
    log.info("Loading Uniprot metadata")
    dbconn = duckdb.connect(args.db, read_only=True)
    uniprot_tbl = (
        dbconn.execute("""
            SELECT rowid, primary_accession, dataset, sequence, names, protein_names, proteinexistence
            FROM uniprot
        """)
        .df()
        .set_index("rowid")
    )
    assert uniprot_tbl['primary_accession'].value_counts().max() == 1, "primary_accession should be a unique key"

    # NB: sequences may be redundant!!
    seq2rowids = uniprot_tbl.index.groupby(uniprot_tbl['sequence'])

    log.info("Building sequence lookup tables from Uniprot sequences")
    # TODO: lookup approach for peptide -> proteins
    kmer_hashes = defaultdict(lambda: set())
    for seq in tqdm(
        seq2rowids.keys(),
        total=len(seq2rowids),
        desc="Building sequence lookup tables",
        unit="sequence",
    ):
        for kmer in window_overlapping(seq, 6, 1):
            kmer_hashes[kmer].update(seq2rowids[seq])

    log.info("Loading input files")
    # iterate through files, read peptides (plus modifications? arent used...)
    peptides_seen = set()
    for fn in tqdm(fns, desc="Reading files"):
        for chunk_tbl in pd.read_table(
            fn, sep="\t", usecols=["Sequence", "Reverse", "Potential contaminant"], chunksize=10_000, low_memory=False,
        ):
            chunk_tbl = chunk_tbl.loc[lambda df: df[["Reverse", "Potential contaminant"]].isna().all(axis=1)]
            peptides_seen.update(chunk_tbl["Sequence"].unique())

    log.info("Got %d distinct peptide sequences" % len(peptides_seen))

    log.info("Building peptide-to-protein graph")
    g = nx.DiGraph()
    for peptide in tqdm(
        peptides_seen, desc="Matching peptides to protein sequences", unit="peptide"
    ):
        candidates = reduce(
            set.intersection,
            (kmer_hashes[kmer] for kmer in window_overlapping(peptide, 6, 1)),
        )

        results = (
            uniprot_tbl.loc[list(candidates)]["sequence"]
            .str.contains(peptide)
            .loc[lambda x: x]
            .index
        )
        for rowid in results:
            # DEBUG
            if peptide == "EHALLAYTQKR":
                print(peptide, candidates, rowid)
            g.add_edge(rowid, peptide)

    rowid2uniprotcost = (
        uniprot_tbl["dataset"].replace({"Swiss-Prot": 0, "TrEMBL": 100}).astype(float)
        + uniprot_tbl["proteinexistence"]
        .replace(
            {
                "evidence at protein level": 0,
                "evidence at transcript level": 5,
                "inferred from homology": 100,
                "predicted": 100,
                "uncertain": 200,
            }
        )
        .astype(float)
    ).to_dict()

    def solve_maxflow(sg, top, bottom) -> dict:
        for edge in sg.edges:
            sg.edges[edge]["capacity"] = 1
            sg.edges[edge]["cost"] = 0

        sg.add_node("Bottom")
        for node in bottom:
            sg.add_edge(node, "Bottom", capacity=1, cost=0)

        sg.add_node("Top")
        for node in top:
            sg.add_edge("Top", node, capacity=100_000, cost=1 + rowid2uniprotcost[node])

        d = nx.algorithms.max_flow_min_cost(
            sg, "Top", "Bottom", capacity="capacity", weight="cost"
        )

        return {k for k in d["Top"] if d["Top"][k] > 0}

    log.info("Assigning peptides to proteins")
    peptide_assignments = {}
    for cc in tqdm(
        nx.connected_components(g.to_undirected()),
        desc="Assigning peptides to proteins",
        unit="Peptide",
    ):
        sg = nx.subgraph(g, cc)
        assert nx.is_bipartite(sg)
        top = [node for node in cc if type(node) is int]
        bottom = [node for node in cc if type(node) is str]
        if any(n == "EHALLAYTQKR" for n in bottom):
            print(cc)

        nr_nodes = solve_maxflow(sg.copy(), top, bottom)

        # Note that the algorithm is non-deterministic and might select a subset of proteins in case of ties
        # Take each other protein that (1) has the same set of peptides assigned and (2) has at least the same score (neg. cost)
        # map set-of-proteins to set-of-peptides
        Z = defaultdict(lambda: set())
        for node in bottom:
            k = frozenset(sg.predecessors(node))
            Z[k].add(node)

        # map set-of-proteins to (set of non-redunant proteins, set of semi-non-redundant proteins)
        M = {}
        for k in Z:
            kp = frozenset(node for node in k if node in nr_nodes)
            m = min(rowid2uniprotcost[node] for node in kp)
            ks = frozenset(
                node for node in k if node in nr_nodes or rowid2uniprotcost[node] <= m
            )
            M[k] = (kp, ks)

        for k in Z:
            kp, ks = M[k]
            for peptide in Z[k]:
                peptide_assignments[peptide] = (k, kp, ks)

    rowid2accession = uniprot_tbl["primary_accession"].to_dict()
    accession2hgncids = (
        dbconn.execute("SELECT primary_accession, hgnc_refs FROM uniprot")
        .df()
        .set_index("primary_accession")['hgnc_refs']
        .loc[lambda x: x.map(len) > 0]
        .map(set)
    ).to_dict()

    peptide_df = pd.DataFrame.from_dict(peptide_assignments, orient='index')
    peptide_uniprot_df = peptide_df.map(
        lambda x: set(rowid2accession[rowid] for rowid in x)
    )
    peptide_hgnc_df = peptide_uniprot_df.map(
        lambda x: reduce(
            set.union,
            (accession2hgncids.get(acc, set()) for acc in x)
        )
    )

    # format and write to file
    (
        peptide_uniprot_df
        .map(lambda x: sorted(x, key=natural_sort_key))
        .map(lambda x: ";".join(x))
        .set_axis(["uniprot_ids", "non_redundant_uniprot_ids", "semi_non_redundant_uniprot_ids"], axis=1)
        .rename_axis("peptide", axis=0)
        .reset_index()
    ).to_csv(outbase / "peptide_assignments.uniprot.tsv.gz", sep="\t", compression="gzip", header=True, index=False)

    (
        peptide_hgnc_df
        .map(lambda x: sorted(x, key=natural_sort_key))
        .map(lambda x: ";".join(x))
        .set_axis(["hgnc_ids", "non_redundant_hgnc_ids", "semi_non_redundant_hgnc_ids"], axis=1)
        .rename_axis("peptide", axis=0)
        .reset_index()
    ).to_csv(outbase / "peptide_assignments.hgnc.tsv.gz", sep="\t", compression="gzip", header=True, index=False)


if __name__ == "__main__":
    main()
