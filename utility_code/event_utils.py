import numpy as np
import pandas as pd


def accept_event(
    *,
    plotdf: pd.DataFrame,
    event_key: str,
    event_start: int,
    event_end: int,
    cn_cutoff: float,
    pt_cutoff: float,
    expected_cn: float,
    **kwargs,
) -> bool:
    """expected_cn can be either arm level information (e.g. median copy number)
    or genome level information (e.g. ploidy)"""
    if isinstance(plotdf, pd.Series):
        plotdf = plotdf.to_frame().T

    # Segments
    seg_start = plotdf["Start.bp"]
    seg_end = plotdf["End.bp"]

    # Add column for tracking nucleotide counts
    plotdf.loc[:, "nt_inside"] = 0

    # CASE: The start of the segment is outside the event boundary but the end is not
    case1 = (seg_end > event_start).values & (seg_start < event_start).values
    # Log the number of nucleotides inside the event boundary
    plotdf.loc[case1, "nt_inside"] = plotdf[case1]["End.bp"].subtract(event_start)

    # CASE: The start and end are contained in the event boundary
    case2 = (seg_end < event_end).values & (seg_start > event_start).values
    # Log the number of nucleotides inside the event boundary
    plotdf.loc[case2, "nt_inside"] = plotdf[case2]["length"]

    # CASE: The start is contained but the end is not
    case3 = (
        (seg_end > event_end).values
        & (seg_start > event_start).values
        & (seg_start < event_end).values
    )
    plotdf.loc[case3, "nt_inside"] = (
        plotdf[case3]["Start.bp"].subtract(event_end).map(abs)
    )

    # Decision
    denominator = plotdf["nt_inside"].sum()
    if event_key == "amp":
        cn_mask = (
            plotdf["expected_total_cn"] > expected_cn + cn_cutoff
        ).values  # meets CN threshold
    else:
        cn_mask = (
            plotdf["expected_total_cn"] < expected_cn - cn_cutoff
        ).values  # meets CN threshold
    nz_mask = (plotdf["nt_inside"] != 0).values  # is non-zero
    numerator = plotdf[cn_mask & nz_mask]["nt_inside"].sum()
    if not denominator:  # avoid division by zero
        return False, 0.0
    # Return the boolean and the fraction itself
    event_fraction = numerator / denominator
    return event_fraction > pt_cutoff, event_fraction


def assign_segments_to_arms(
    *, cdat: pd.DataFrame, centromeres: pd.DataFrame, chromosome: int
) -> pd.DataFrame:
    """Updates cdat with chromosome arm labels and truncates segments to one or the other arm

    Split / Truncate / Drop Segments:
        A segment may be only on the p arm, only on the q arm, span both arms, or
        be inside the centromere entirely. In each case, we truncate either end that
        is inside the centromere and ignore entirely centromere exclusive segments.
    """
    # Centromere start and end locations
    C_start = centromeres.loc[chromosome]["start"]
    C_end = centromeres.loc[chromosome]["end"]
    # Segment start and end locations
    seg_start = cdat["Start.bp"]
    seg_end = cdat["End.bp"]

    tmp = cdat.copy().reset_index()
    # p arm only
    case1 = (
        (seg_start < C_start).values
        & (seg_end > C_start).values
        & (seg_end <= C_end).values
    )
    tmp.loc[case1, "End.bp"] = C_start  # clip end
    seg_end = tmp["End.bp"]

    case2 = (seg_start < C_start).values & (seg_end <= C_start).values
    tmp.loc[case2, "Chromosome"] = f"{chromosome}p"  # assign to p arm

    # q arm only
    case3 = (seg_start > C_start).values & (
        seg_start <= C_end
    ).values  # inside centromere
    tmp.loc[case3, "Start.bp"] = C_end  # clip start
    seg_start = tmp["Start.bp"]

    case4 = (seg_start >= C_end).values  # set to q arm
    tmp.loc[case4, "Chromosome"] = f"{chromosome}q"

    # null case -> segment entirely inside centromere
    case5 = (
        (seg_start > C_start)
        & (seg_start < C_end)
        & (seg_end > C_start)
        & (seg_end < C_end)
    ).values  # drop
    if not tmp[case5].empty:
        tmp = tmp.drop(tmp[case5].index)
        seg_start = tmp["Start.bp"]
        seg_end = tmp["End.bp"]

    # spanning case - split segment into p arm and q arm
    case6 = (seg_start < C_start).values & (
        seg_end > C_end
    ).values  # split segments and clip
    split = tmp[case6]
    if not split.empty:
        if "index" in split.columns:
            split = split.drop("index", axis=1)
            tmp = tmp.drop("index", axis=1)

        if len(split) > 1:
            p = pd.DataFrame(
                [
                    [f"{chromosome}p"] * len(split),  # chromosome arm
                    split["Start.bp"].tolist(),  # start bp
                    [C_start] * len(split),  # end bp
                    [0] * len(split),  # length placeholder
                    split["expected_total_cn"].tolist(),  # exp cn
                ]
            ).T
            q = pd.DataFrame(
                [
                    [f"{chromosome}q"] * len(split),  # chromosome arm
                    [C_end] * len(split),  # start bp
                    split["End.bp"].tolist(),  # end bp
                    [0] * len(split),  # length placeholder
                    split["expected_total_cn"].tolist(),  # exp cn
                ]
            ).T
            pq = pd.concat([p, q])
            pq.columns = tmp.columns
            pq.loc[:, "length"] = pq["End.bp"].subtract(pq["Start.bp"])
            pq.index = split.index.tolist() * 2
        else:
            p = [
                f"{chromosome}p",  # chromosome arm
                split["Start.bp"].tolist()[0],  # start bp
                C_start,  # end bp
                C_start - split["Start.bp"].tolist()[0],  # length placeholder
                split["expected_total_cn"].tolist()[0],  # exp cn
            ]

            q = [
                f"{chromosome}q",  # chromosome arm
                C_end,  # start bp
                split["End.bp"].tolist()[0],  # end bp
                split["End.bp"].tolist()[0] - C_end,  # length placeholder
                split["expected_total_cn"].tolist()[0],  # exp cn
            ]
            pq = pd.concat([pd.Series(p), pd.Series(q)], axis=1).T
            pq.columns = tmp.columns
            pq.index = split.index.tolist() * len(pq)
        tmp = (
            pd.concat([tmp.drop(split.index), pq])
            .sort_values("Start.bp")
            .set_index("Chromosome")
        )
    if tmp.index.name != "Chromosome":
        tmp = tmp.set_index("Chromosome")
    if "index" in tmp.columns:
        tmp = tmp.drop("index", axis=1)
    return tmp
