import pandas as pd
import duckdb


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("-d", "--db", metavar="OMNI_GENE_DUCKDB", required=True)
    ap.add_argument("-p", "--peptideinfo", required=True)
    ap.add_argument("-m", "--metadata", required=True)
    ap.add_argument("infile", nargs=1)

    args = ap.parse_args()

    dbconn = duckdb.connect(args.db, read_only=True)

    hgnc_tbl = dbconn.execute("SELECT hgnc_id, symbol, name, locus_group, locus_type, ensembl_gene_id FROM hgnc").df()
    dd = hgnc_tbl.set_index("hgnc_id").to_dict(orient="index")

    peptide_assignments = pd.read_csv(args.infile[0], sep="\t", low_memory=False)
    hgnc2pep = (
        peptide_assignments
        .fillna({'semi_non_redundant_hgnc_ids': "N/A"})
        .groupby('semi_non_redundant_hgnc_ids')
        ['peptide']
        .agg(set)
    )
    peptide_info = pd.DataFrame({
        "peptide_sequences": hgnc2pep.map(lambda s: ";".join(sorted(s))),
        "peptide_count": hgnc2pep.map(len),
    })
    peptide_info.to_csv(args.peptideinfo, sep="\t", index=True, header=True)

    # join columns of data with semicolons, create new rows aligning to semicolon-joined HGNC IDs:
    output = pd.DataFrame.from_dict(
        {
            k: {
                col + "s": ";".join(dd[hgncid][col] for hgncid in k.split(";"))
                for col in hgnc_tbl.iloc[:, 1:].columns
            }
            for k in set(hgnc2pep.index) - {"N/A"}
        }, orient='index'
    )

    output.rename_axis(index='hgnc_ids').reset_index().sort_values("hgnc_ids").to_csv(args.metadata, sep="\t", index=False, header=True)


if __name__ == "__main__":
    main()
