#!/usr/bin/env python
# coding: utf-8
import numpy as np
import pandas as pd
import xarray
import matplotlib
from matplotlib import pyplot as plt
matplotlib.use("pdf")
import seaborn

np.random.seed(42)

BASE_YLIM = (-2, 3)
BASE_COEFLIM = (-0.25, 0.25)
MISSING_COLOR = 'darkred'


sample_subset = pd.read_csv("pipeline_activation/samplesheet.tsv", sep="\t")

intensities = pd.read_csv(
    "pipeline_activation/intensity-values.tsv",
    sep="\t", index_col=0,
)

w_valid = intensities['w_valid']

Y = intensities[[col for col in sample_subset['sample_label']]].values
with np.errstate(invalid='ignore', divide='ignore'):
    Yl = np.log10(np.ma.masked_less_equal(Y, 0))

log_intensity_medians = pd.read_csv("pipeline_activation/model-output/log_intensity_medians.tsv", sep="\t", index_col=0)
log_intensity_medians = log_intensity_medians.values.squeeze()

Yln = Yl - (log_intensity_medians / np.log(10))

assert np.isclose([
    np.quantile(Yln[w_valid][:, i].compressed(), q=0.5, axis=0)
    for i in range(Yln.shape[1])
], 0.).all()


def round_away(x, decimals=0):
    f = 10 ** decimals
    return np.copysign(np.ceil(np.abs(x) * f) / f, x)


lim = np.abs(Yln).max()

map_ds = xarray.load_dataset("pipeline_activation/model-output/global_map.h5")
E = xarray.concat(
    [map_ds['E'], map_ds['Ederiv'].rename({"deriv_coef": "coef"})],
    dim='coef',
)
hdi_ds = xarray.open_dataset(
    "pipeline_activation/model-output/pergene_hdis.h5"
)
assert 'H90' in list(hdi_ds.variables)
assert 'Hd90' in list(hdi_ds.variables)
H = xarray.concat(
    [hdi_ds['H90'], hdi_ds['Hd90'].rename({"deriv_coef": "coef"})],
    dim='coef',
)
assert w_valid.sum() == E.shape[0]
assert (E.coords["coef"] == H.coords["coef"]).all()


plot_coefs = [
    'activation',
    'RBPness',
    'act_deltaRBPness',
    'activation_FP',
]
assert all(coef in E.coords['coef'] for coef in plot_coefs)


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("-m", "--metadata", required=True)
    ap.add_argument("-o", "--outfile", required=True)
    ap.add_argument("--only-libtype", type=str, choices=['OOPS', 'fullproteome'], required=False)
    ap.add_argument("--with-coef", action="store_true")
    ap.add_argument("--with-legend", action="store_true")
    ap.add_argument("query")
    args = ap.parse_args()

    hgnc_metadata = pd.read_csv(args.metadata, sep="\t", header=0, index_col=0)
    assert (hgnc_metadata.index == intensities.index).all()

    # First, try exact match
    has_symbol = (
        hgnc_metadata['symbols'] == args.query
    )
    if has_symbol.sum() == 1:
        assert has_symbol[w_valid].sum() == 1
    else:
        has_symbol = (
            hgnc_metadata['symbols']
            .str.split(";").apply(set)
            .apply(lambda s: args.query in s)
        )

        assert has_symbol.sum() == 1
        assert has_symbol[w_valid].sum() == 1

    iidx = has_symbol[w_valid].argmax()

    # Subset samples if necessary
    if args.only_libtype is not None:
        sample_subset = sample_subset.groupby(['library_type']).get_group((args.only_libtype, ))
        sample_idxs = sample_subset.index
        d = Yln[w_valid][:, sample_idxs][iidx]
    else:
        d = Yln[w_valid][iidx]

    # Determine limits
    if d.min() < BASE_YLIM[0]:
        print("WARNING: changing raw data ymin")
        YMIN = round_away(d.min(), 1)
    else:
        YMIN = BASE_YLIM[0]
    if d.max() > BASE_YLIM[1]:
        print("WARNING: changing raw data ymax")
        YMAX = round_away(d.max(), 1)
    else:
        YMAX = BASE_YLIM[1]

    # Setup marker/plotting data
    MARKERS = dict(zip(
        sorted(sample_subset["donor_id"].unique()),
        "osD",
    ))

    marker_data = {}
    for igroup, (groupkey, g) in enumerate(sample_subset.groupby(['library_type', 'activated', 'crosslinked'])):
        for irow, idx in enumerate(g.sort_values("donor_id").index):
            color_idx = (
                0 if groupkey[2] == "yes"
                else 1
            ) + (
                6 if groupkey[1] == "yes"
                else 0
            )
            color = seaborn.color_palette("tab20")[color_idx]
            if groupkey[0] == "OOPS":
                linecolor = color
                fillcolor = color
            else:
                linecolor = color
                fillcolor = "white"

            marker_data[idx] = {
                "xpos": 0.5 * (igroup + 0.8 * (-0.5 + (irow / (len(g) - 1)))),
                "linecolor": linecolor,
                "fillcolor": fillcolor,
                "marker": MARKERS[g.loc[idx]['donor_id']],
            }

    if args.with_coef:
        v = H.isel(gene=iidx).sel(coef=plot_coefs)

        if v.min().item() < BASE_COEFLIM[0]:
            print("WARNING: changing coef ymin")
            COEF_YMIN = round_away(v.min().item(), 1)
        else:
            COEF_YMIN = BASE_COEFLIM[0]
        if v.max().item() > BASE_COEFLIM[1]:
            print("WARNING: changing coef ymax")
            COEF_YMAX = round_away(v.max().item(), 1)
        else:
            COEF_YMAX = BASE_COEFLIM[1]

    # Generate the plot
    ncols = (2 if args.with_coef else 1)
    figwidth = (7 if args.with_coef else 4.6)

    fig, axs = plt.subplots(
        1, ncols, figsize=(figwidth, 6),
        dpi=100,
        gridspec_kw=(
            dict(width_ratios=[2, 1])
            if args.with_coef
            else {}
        ),
        sharex=False, sharey=False, squeeze=False,
        layout='none',
    )
    fig.subplots_adjust(bottom=0.25, wspace=0.5)
    fig.set_layout_engine('none')
    if args.with_coef:
        axdata, axcoef = axs[0]
    else:
        axdata = axs[0][0]

    # Scatter points
    scs = []
    for datum, row in zip(d.filled(YMIN), pd.DataFrame.from_dict(marker_data, orient='index').sort_index().itertuples()):
        row = row._asdict()
        sc = axdata.scatter(
            [row['xpos']], [datum],
            marker=row['marker'],
            s=5 ** 2,
            zorder=10.,
            lw=2, color=row['fillcolor'], edgecolors=row['linecolor'],
        )
        if datum <= YMIN:
            sc.set_facecolor(MISSING_COLOR)
        sc.set_clip_on(False)
        scs.append(sc)

    # x-axis labeling:
    axdata.set_xticks(np.arange((4 if args.only_libtype is None else 2)) + 0.5 / 2)
    axdata.set_xticklabels([
        "%s\n%s" % (x[0], {"yes": "activated", "no": "resting"}[x[1]])
        for x in
        list(sample_subset.groupby(['library_type', 'activated']).groups.keys())
    ], rotation=90,
    )

    axdata.set_ylabel("Centered log-intensity")
    axdata.set_title(args.query)

    axdata.set_ylim(YMIN, 3)
    seaborn.despine(ax=axdata)

    # y-axis and missing value indicator:
    hl = axdata.axhline(YMIN, color=MISSING_COLOR, zorder=9, lw=2.)
    hl.set_clip_on(False)
    axdata.yaxis.set_major_locator(plt.matplotlib.ticker.MultipleLocator(1))
    axdata.yaxis.set_major_formatter(plt.matplotlib.ticker.FuncFormatter(lambda x, pos: (x if x > YMIN else "ND")))

    if args.with_coef:
        axcoef.errorbar(
            np.arange(len(plot_coefs)),
            y=E.isel(gene=iidx).sel(coef=plot_coefs).values,
            yerr=[
                E.isel(gene=iidx).sel(coef=plot_coefs).values - H.isel(gene=iidx).sel(coef=plot_coefs).values[:, 0],
                H.isel(gene=iidx).sel(coef=plot_coefs).values[:, 1] - E.isel(gene=iidx).sel(coef=plot_coefs).values,
            ],
            fmt='none',
            color='black',
            lw=3., elinewidth=3., capthick=3., capsize=3.,
            zorder=10.,
        )
        axcoef.scatter(
            x=np.arange(len(plot_coefs)),
            y=E.isel(gene=iidx).sel(coef=plot_coefs).values,
            zorder=11.,
            edgecolor='black', color='white', lw=2.,
        )
        axcoef.set_xlim(-0.5, 0.5 + len(plot_coefs) - 1)
        axcoef.set_xticks(np.arange(len(plot_coefs)))
        axcoef.set_xticklabels(plot_coefs, rotation=90)
        axcoef.axhline(0., lw=2., ls="--", color="black", zorder=2)

        axcoef.set_ylim(COEF_YMIN, COEF_YMAX)

        axcoef.set_ylabel("Coefficient\n(90% credible interval)")
        axcoef.set_title(args.query)

    # # slightly rotate x-axis labels
    # for ax in (axdata, axcoef):
    #     plt.setp(ax.get_xticklabels(), rotation=60, ha="right", va="top")

    # Legend box:
    if args.with_legend:
        axdata.scatter([0.05], [0.95], s=5 ** 2, marker='o', color='black', transform=axdata.transAxes, gid="legend")
        axdata.scatter([0.09], [0.95], s=5 ** 2, marker='o', color='grey', transform=axdata.transAxes, gid="legend")
        axdata.text(0.16, 0.95, "CL / nCL", ha="left", va="center", transform=axdata.transAxes, gid="legend")

        axdata.scatter([0.05], [0.90], s=5 ** 2, marker='o', color=seaborn.color_palette("tab20")[6], transform=axdata.transAxes, gid="legend")
        axdata.scatter([0.09], [0.90], s=5 ** 2, marker='o', color=seaborn.color_palette("tab20")[0], transform=axdata.transAxes, gid="legend")
        axdata.text(0.16, 0.90, "Act / Rest", ha="left", va="center", transform=axdata.transAxes, gid="legend")

        axdata.scatter([0.05], [0.85], s=5 ** 2, marker='o', color='black', transform=axdata.transAxes, gid="legend")
        axdata.scatter([0.09], [0.85], s=5 ** 2, marker='o', color='white', edgecolor='black', transform=axdata.transAxes, gid="legend")
        axdata.text(0.16, 0.85, "OOPS / TCL", ha="left", va="center", transform=axdata.transAxes, gid="legend")

        axdata.scatter([0.05], [0.80], s=5 ** 2, marker='o', color='black', transform=axdata.transAxes, gid="legend")
        axdata.scatter([0.09], [0.80], s=5 ** 2, marker='s', color='black', transform=axdata.transAxes, gid="legend")
        axdata.scatter([0.13], [0.80], s=5 ** 2, marker='D', color='black', transform=axdata.transAxes, gid="legend")
        axdata.text(0.16, 0.80, "donor pool", ha="left", va="center", transform=axdata.transAxes, gid="legend")

        axdata.add_patch(
            plt.matplotlib.patches.Rectangle((0., 1.,), 0.45, -0.25, transform=axdata.transAxes, facecolor='none', edgecolor='black', lw=2., gid="legend"),
        )

    fig.savefig(args.outfile, format="pdf")
