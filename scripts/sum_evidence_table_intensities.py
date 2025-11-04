import sys
import pandas as pd

peptide_assignments = pd.read_csv(sys.argv[1], sep="\t", low_memory=False)
assert peptide_assignments['peptide'].value_counts().max() == 1

dtbl = pd.read_csv(sys.argv[2], sep="\t", low_memory=False)

lut = peptide_assignments.set_index("peptide").iloc[:, 2].fillna('N/A').to_dict()
gene_coords = sorted(set(lut.values()))
wtbl = (
    dtbl['Intensity']
    .fillna(0.)
    .groupby([dtbl['Sequence'].map(lut.get).rename("Gene"), dtbl['Experiment']])
    .agg("sum")
    .to_frame().reset_index()
    .pivot(columns='Experiment', index='Gene', values='Intensity')
    .reindex(gene_coords)
)

wtbl.to_csv(sys.stdout, sep="\t", index=True, header=True)
