import os
import numpy as np
import pandas as pd
import dalmatian
from .event_utils import assign_segments_to_arms, accept_event


def truth_table(configs: list) -> pd.DataFrame:
    tt = list()
    for i in range(len(configs)):
        arm = configs[i]["arm"]
        plotdfs = configs[i]["plotdfs"]
        plotdfs = plotdfs[~plotdfs.index.duplicated(keep="first")]  # dedup index
        plotdfs = plotdfs[["accepted"]]
        plotdfs.columns = [f"{a}_{arm}" for a in plotdfs.columns]
        tt.append(plotdfs)
    tt = pd.concat(tt, axis=1)
    return tt


def ploidy_mtx(
    namespace: str, workspace: str, sample_event_data: pd.DataFrame
) -> pd.Series:
    P = dalmatian.WorkspaceManager(f"{namespace}/{workspace}").get_pairs()
    P = P[~P["absolute_called_ploidy"].isna()]
    samples = sample_event_data.index.get_level_values(0).unique()
    ploidy_mtx = P.loc[[f"{c}_pair" for c in samples]]["absolute_called_ploidy"]
    ploidy_mtx.index = ploidy_mtx.index.str.replace("_pair", "")
    ploidy_mtx = ploidy_mtx.astype(float)
    return ploidy_mtx


def event_data(events: list, lesions: pd.DataFrame) -> pd.DataFrame:
    key = "wide peak boundaries"
    event_data = lesions.set_index("cytoband").T[[key]].copy()
    if any(event_data.index.str.contains("Unnamed")):
        event_data = event_data.drop(
            event_data.index[event_data.index.str.contains("Unnamed")]
        )

    # Add chromosome information
    event_data.loc[:, "chrom"] = (
        event_data[key].str.split(":").str[0].str.replace("chr", "")
    )

    if 'X' in list(event_data['chrom']): #removing X chromosome events
        event_data.drop(event_data[event_data['chrom']=='X'].index, inplace=True)

    event_data["chrom"] = event_data["chrom"].astype(int)

    # Add arm information
    event_data.loc[:, "arm"] = event_data["chrom"].astype(str) + pd.Series(
        event_data.index.str.extract("(q|p)")[0].tolist(), index=event_data.index
    )
    # Add event start boundary
    event_data.loc[:, "start"] = (
        event_data[key].str.split(":").str[-1].str.split("-").str[0]
    )
    event_data["start"] = event_data["start"].astype(int)

    # Add event end boundary
    event_data.loc[:, "end"] = (
        event_data[key].str.split(":").str[-1].str.split("-").str[1]
    )
    event_data["end"] = event_data["end"].astype(int)
    event_data = event_data.sort_index()

    return event_data


def cytoref(cytoband_file: str):
    # Create modified cytoband ref file that includes the indiv. cytoband names
    cytoref = pd.read_csv(cytoband_file, sep="\t")
    cytoref = cytoref.T.reset_index().T  # fix formatting
    cytoref.index = [i for i in range(len(cytoref))]
    cytoref.columns = ["chrom", "start", "end", "band", "stain"]
    cytoref.loc[:, "cytoband"] = (
        cytoref["chrom"].str.replace("chr", "") + cytoref["band"]
    )
    cytoref = cytoref.set_index("cytoband")
    return cytoref


def centromere_coords(cytoband: pd.DataFrame) -> pd.DataFrame:
    centromeres = pd.unique(
        [int(c.replace("p", "").replace("q", "")) for c in cytoband.index]
    )
    centromeres = pd.DataFrame(index=centromeres)
    centromeres.loc[:, "start"] = None
    centromeres.loc[:, "end"] = None
    for c in centromeres.index:
        coords = cytoband.loc[[f"{c}p", f"{c}q"]]
        centromeres.loc[c, "start"] = coords["end"][0]
        centromeres.loc[c, "end"] = coords["start"][-1]
    return centromeres


def sample_event_data(
    sub_segs: pd.DataFrame, centromeres: pd.DataFrame
) -> pd.DataFrame:
    sample_event_data = dict()
    samples = sub_segs.index.get_level_values(0).unique().tolist()
    for sample in samples:
        dats = dict()
        dat = sub_segs.loc[sample].set_index("Chromosome")
        for chromosome in dat.index.unique():
            dats[chromosome] = dict()
            cdat = dat.loc[chromosome]  # chromosome level information
            if isinstance(cdat, pd.Series):
                cdat = (
                    cdat.to_frame()
                    .T.reset_index()
                    .rename(columns={"index": "Chromosome"})
                )
            cdat = assign_segments_to_arms(
                cdat=cdat, centromeres=centromeres, chromosome=chromosome
            )
            dats[chromosome] = cdat
        sample_event_data[sample] = dats
    sample_event_data = (
        pd.concat({k: pd.concat(v) for k, v in sample_event_data.items()})
        .reset_index()
        .drop("level_1", axis=1)
        .set_index(["level_0", "Chromosome"])
    )
    sample_event_data.index.names = ["sample", "Chromosome"]
    return sample_event_data


def copy_number_per_arm(
    sample_event_data: pd.DataFrame, event_data: pd.DataFrame, method: str
) -> pd.DataFrame:
    samples = sample_event_data.index.get_level_values(0).unique().tolist()
    exp_cn_sample_arm = {sample: {} for sample in samples}
    for sample in samples:
        for arm in event_data["arm"]: #1p, 1q, 2p ...
            dat = sample_event_data.loc[sample]
            if arm not in dat.index:
                continue
            dat = dat.loc[arm]
            if isinstance(dat, pd.DataFrame):
                if method == "mean":
                    exp_cn_arm = (
                        dat["length"]
                        .divide(dat["length"].sum())
                        .mul(dat["expected_total_cn"])
                        .sum()
                    )
                elif method == "median":
                    dat = dat.sort_values("expected_total_cn").copy()
                    dat.loc[:, "cumulative_length"] = dat["length"].cumsum()
                    median_length = dat["cumulative_length"].iloc[-1] // 2
                    idx = np.searchsorted(dat["cumulative_length"], median_length)
                    exp_cn_arm = dat.iloc[idx]["expected_total_cn"]
            else:  # single item -> series
                exp_cn_arm = int(round(float(dat["expected_total_cn"]), 1))
            exp_cn_sample_arm[sample][arm] = exp_cn_arm
    # NaN indicates arm was skipped
    exp_cn_sample_arm = pd.DataFrame.from_dict(exp_cn_sample_arm)
    exp_cn_sample_arm = exp_cn_sample_arm[sorted(exp_cn_sample_arm.columns)]
    return exp_cn_sample_arm


def sub_seg_file(absolute_files: str) -> pd.DataFrame:
    """Generate a subset of concatenated absolute seg files"""
    cols = ["Chromosome", "Start.bp", "End.bp", "length", "expected_total_cn", "LOH"]
    files = list()
    samples = list()
    for i, file in enumerate(os.listdir(absolute_files)):
        sample = file.split(".")[0].strip("_pair")
        samples.append(sample)
        files.append(os.path.join(absolute_files, file))

    segs = pd.concat(
        {
            sample: pd.read_csv(file, sep="\t").drop("sample", axis=1)
            for sample, file in zip(samples, files)
        }
    )
    # Make subset sample level seg files
    sub_segs = segs[cols].copy()
    return sub_segs


def plot_configs_focal_level(
    *,
    sample_event_data: pd.DataFrame,
    exp_cn_sample_arm: pd.DataFrame,
    event_data: pd.DataFrame,
    cytoref: pd.DataFrame,
    cytoband: pd.DataFrame,
    cn_cutoff: float,
    pt_cutoff: float,
    event_key: str,
    method: str,
    out_dir: str,
) -> pd.DataFrame:
    configs = list()
    event_fractions = dict()
    # Use the event-data to generate grid plots for each event on each chromosome
    for event, data in event_data.iterrows():

        event_fractions[event] = dict()
        event_start, event_end = data["start"], data["end"]
        chrom, arm = data["chrom"], data["arm"]

        plotdfs = {}
        # Subset to just those samples that had an event
        event_samples = sample_event_data.index.get_level_values(0).unique()
        for sample in event_samples:
            event_fractions[event][sample] = None
            plotdf = sample_event_data.loc[sample].loc[arm]
            if isinstance(plotdf, pd.Series):
                plotdf = plotdf.to_frame().T
            expected_cn = exp_cn_sample_arm[sample].loc[arm]
            accepted, event_fraction = accept_event(
                plotdf=plotdf,
                event=event,
                event_key=event_key,
                event_start=event_start,
                event_end=event_end,
                cn_cutoff=cn_cutoff,
                pt_cutoff=pt_cutoff,
                expected_cn=expected_cn,
                )
            plotdf.loc[:, "accepted"] = accepted
            plotdf.loc[:, "event_fraction"] = event_fraction
            event_fractions[event][sample] = event_fraction
            plotdfs[sample] = plotdf
        plotdfs = pd.concat(plotdfs).sort_index()
        plotdfs.index = plotdfs.index.get_level_values(0)
        configs.append(
            dict(
                plotdfs=plotdfs,
                exp_cn_sample_arm=exp_cn_sample_arm,
                cytoband=cytoband,
                chrom=chrom,
                arm=arm,
                event_start=event_start,
                event_end=event_end,
                event=event,
                out_dir=out_dir,
                method=method,
                cn_cutoff=cn_cutoff,
                pt_cutoff=pt_cutoff,
                event_fractions=event_fractions,
            )
        )
    return configs


def event_fractions(
    saved_configs: dict, cn_cutoffs: list, pt_cutoffs: list
) -> pd.DataFrame:
    s = dict()
    for cn in cn_cutoffs:
        s[cn] = dict()
        for pt in pt_cutoffs:
            df = (
                pd.DataFrame.from_dict(saved_configs[cn][pt], orient="index")
                .fillna(0)
                .T
            )

            t = df[[c for c in df.columns if "accepted" not in c]] > pt
            [
                df.insert(
                    i,
                    f"{df.iloc[:, :i].columns[-1]}_accepted",
                    t[df.iloc[:, :i].columns[-1]].values.flatten(),
                )
                for i in range(df.shape[1], 0, -1)
            ]  # inplace truth column insertion
            s[cn][pt] = df
    s = pd.concat({k: pd.concat(v) for k, v in s.items()})
    return s


def acceptance_rates(event_fractions: pd.DataFrame) -> pd.DataFrame:
    cn_cutoffs = event_fractions.index.get_level_values(0).unique().tolist()
    pt_cutoffs = event_fractions.index.get_level_values(1).unique().tolist()
    sums = dict()
    for cn in cn_cutoffs:
        sums[cn] = dict()
        for pt in pt_cutoffs:
            data = event_fractions.loc[cn].loc[pt]
            data = data[data.columns[data.columns.str.contains("accepted")]]
            n_samples = len(data.index)
            data = data.sum()
            data.name = pt
            data = data.divide(n_samples)
            sums[cn][pt] = data
    sums = pd.concat(
        {k: pd.DataFrame.from_dict(v, orient="index") for k, v in sums.items()}
    )
    sums.columns = sums.columns.str.replace("_accepted", "")
    return sums


def event_list_to_dataframe(event_list, cytoref):
    '''
    takes list of events and returns a dataframe of individual events
    '''
    event_list = pd.Series(event_list)

    for i in range(len(event_list)):
        if event_list[i].find('_') > 0:
            event_list[i] = event_list[i][0:-1]

    df_events = pd.DataFrame({'event':event_list.str.extract("(\w{4})")[0],
                            'cytoband':event_list.str.split("gain|loss").str[1]})

    df_split = df_events['cytoband'].str.split('(p|q)')
    df_events['chrom'] = df_split.str[0]
    df_events['p/q'] = df_split.str[1]
    df_events['arm'] = df_events['chrom'] + df_events['p/q']

    df_single = df_events[~df_events['cytoband'].str.contains('-')]
    locs = df_events['cytoband'].str.split('(p|q)').str[2]
    df_single['start_cyt'], df_single['end_cyt'] = locs, locs

    df_multiple = df_events[df_events['cytoband'].str.contains('-')]
    locations = df_multiple['cytoband'].str.split('(p|q)').str[2]
    starts = locations.str.split('(-)').str[0]
    ends = locations.str.split('(-)').str[2]
    df_multiple['start_cyt'], df_multiple['end_cyt'] = starts, ends

    df_events = pd.concat([df_single, df_multiple], ignore_index=True)

    start_bps, end_bps = [], []
    for i in df_events.index:
        row = df_events.loc[i]
        start_cyt = row['chrom'] + row['p/q'] + row['start_cyt']
        end_cyt = row['chrom'] + row['p/q'] + row['end_cyt']

        start_bps.append(cytoref.loc[start_cyt]['start'])
        end_bps.append(cytoref.loc[start_cyt]['end'])

    df_events['start'] = start_bps
    df_events['end'] = end_bps

    return df_events













###
