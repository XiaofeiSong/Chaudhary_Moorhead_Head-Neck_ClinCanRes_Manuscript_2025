import os
import pandas as pd
import matplotlib.pyplot as plt

plt.style.use("default")


def plot(
    *,
    plotdfs: pd.DataFrame,
    cytoband: pd.DataFrame,
    exp_cn_sample_arm: pd.DataFrame = None,
    chrom: int = None,
    arm: str = None,
    event: str = None,
    event_start: int = None,
    event_end: int = None,
    out_dir: str = None,
    method: str = None,
    cn_cutoff: float = None,
    pt_cutoff: float = None,
    show: bool = False,
    **kwargs,
) -> None:
    """Generates a grid of plots for a given chromosomal event"""
    if event is None and arm is not None:
        event = arm
    # Centromere information
    if chrom is None:
        chrom = arm.replace("p", "").replace("q", "")
    coords = cytoband.loc[[f"{chrom}p", f"{chrom}q"]]
    ## | p_end -> | <- centromere -> | <- q_start |
    centromere = (coords.iloc[0]["end"], coords.iloc[1]["start"])

    # Clip event boundaries to make the plot look neater
    # This doesn't affect event selection (upstream of plotting)
    S, E, C0, C1 = event_start, event_end, centromere[0], centromere[1]
    if S < C0 and E > C1:  # event overlaps the centromere
        if arm == "p":
            E = C0
        else:
            S = C1
    elif S > C0 and S < C1 and E > C1:  # event starts in centromere
        S = C1
    elif S < C0 and S < C1 and E > C0:  # event ends in centromere
        E = C0

    samples = plotdfs.index.unique()
    rows, cols = int(len(samples) / 2) + (len(samples) % 2 > 0), 2
    width = 15
    height = rows * 4
    fig, axes = plt.subplots(rows, cols, figsize=(width, height), sharex=True, sharey=True)
    print(f"Start plotting {event}...")
    for ax, samp in zip(axes.flatten(), samples):
        plotdf = plotdfs.loc[samp]
        if isinstance(plotdf, pd.Series):
            plotdf = plotdf.to_frame().T
        # Plot the individual segments
        for i, (row, dat) in enumerate(plotdf.iterrows()):
            x = (dat["Start.bp"], dat["End.bp"])
            y = (dat["expected_total_cn"], dat["expected_total_cn"])
            if i == len(plotdf) - 1:
                ax.plot(x, y, color="b", linewidth=1, label="Segment")
            else:
                ax.plot(x, y, color="b", linewidth=1)
            ax.set_title(
                f"{samp} | Event: {event} | Accepted: {plotdf['accepted'].unique()[-1]}",
                fontsize=12,
            )
            # Plot the event band boundaries
            if i == len(plotdf) - 1 and event != arm:
                ax.axvspan(
                    S,  # event_start
                    E,  # event_end
                    alpha=0.20,
                    color="k",
                    linestyle="--",
                    hatch="/",
                    label=event,
                )
        # Plot the centromere location
        ax.axvspan(
            C0,
            C1,
            alpha=0.20,
            color="b",
            linestyle="--",
            hatch="/",
            label="Centromere",
        )
        xmin = plotdf["Start.bp"].min() if arm[-1] == "p" else centromere[1]
        xmax = plotdf["End.bp"].max() if arm[-1] == "q" else centromere[0]
        # Plot the ploidy, if present in plotdf
        if "ploidy" in plotdf.columns:
            ploidy = plotdf["ploidy"].unique()[-1]
            ploidy_label = "wxs_ploidy"
            ax.hlines(
                y=ploidy,
                xmin=xmin,
                xmax=xmax,
                ls="--",
                color="green",
                linewidth=2,
                alpha=0.75,
                label=ploidy_label,
            )
        if exp_cn_sample_arm is not None:
            # Plot the expected/median copy number the sample + arm
            ecna = exp_cn_sample_arm[samp].loc[arm]
            ax.hlines(
                y=ecna,
                xmin=xmin,
                xmax=xmax,
                ls="--",
                color="red",
                linewidth=2,
                alpha=0.75,
            )
            ecna_label = r"$E[CN_{arm}]$" if method == "mean" else r"$Median[CN_{arm}]$"
            ax.hlines(
                ecna,
                ls="--",
                xmin=xmin,
                xmax=xmax,
                color="orange",
                linewidth=1.5,
                alpha=0.75,
                label=ecna_label,
            )
    print(f"End plotting {event}...")
    ax.legend()
    out_dir = f"{out_dir}/cn-cut-{cn_cutoff}_pct-cut-{pt_cutoff}"
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    print(f"Saving {event}...")
    fig.savefig(f"{out_dir}/{event}.pdf", bbox_inches="tight")
    print("Done...\n\n")
    plt.show() if show else plt.close()
