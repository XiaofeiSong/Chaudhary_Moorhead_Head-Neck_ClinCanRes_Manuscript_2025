import numpy as np
import pandas as pd
import os
import itertools
# import event_utils
from event_utils import assign_segments_to_arms, accept_event

_chrmap = dict(
    zip(["chr" + str(x) for x in list(range(1, 23)) + ["X", "Y"]], range(1, 25))
)


def parse_cytoband(cytoband: str, annotate_arms: bool = False):
    # some cytoband files have a header, some don't; we need to check
    has_header = False
    with open(cytoband, "r") as f:
        if f.readline().startswith("chr\t"):
            has_header = True

    cband = pd.read_csv(
        cytoband,
        sep="\t",
        names=["chr", "start", "end", "band", "stain"] if not has_header else None,
    )
    cband["chr"] = cband["chr"].apply(lambda x: _chrmap[x])

    chrs = cband["chr"].unique()
    ints = dict(zip(chrs, [{0} for _ in range(0, len(chrs))]))
    last_end = None
    last_stain = None
    last_chrom = None
    for _, chrom, start, end, _, stain in cband.itertuples():
        if start == 0:
            if last_end is not None:
                ints[last_chrom].add(last_end)
        if stain == "acen" and last_stain != "acen":
            ints[chrom].add(start)
        if stain != "acen" and last_stain == "acen":
            ints[chrom].add(start)

        last_end = end
        last_stain = stain
        last_chrom = chrom
    ints[chrom].add(end)

    CI = np.full([len(ints), 4], 0)
    for c in chrs:
        CI[c - 1, :] = sorted(ints[c])

    cytoband = pd.DataFrame(
        np.c_[np.tile(np.c_[np.r_[1:25]], [1, 2]).reshape(-1, 1), CI.reshape(-1, 2)],
        columns=["chr", "start", "end"],
    )
    if annotate_arms: # with p, q information and length
        cytoband.loc[:, "chr"] = (
            cytoband["chr"].astype(str)
            + pd.Series(["p", "q"]).iloc[cytoband.index % 2].values
        ).values
        cytoband.loc[:, "length"] = cytoband["end"].subtract(cytoband["start"])
    cytoband = cytoband.set_index("chr")
    return cytoband



def create_categorical_comut_df(data, data_type):
    df = pd.DataFrame({'sample':data['sample'], 'category':[data_type for q in range(len(data))], 'value':data[data_type]})
    df.drop_duplicates(inplace=True)

    return df




def add_most_recent_pretreatment_to_df(df):
    '''
    assumes sample column as index without name and merged annotation file
    '''

    df_subset = df[(df['indicator'].str.contains('tiss'))|(df['indicator'].str.contains('S'))|(df['indicator'].str.contains('C1D1'))]

    for part in df_subset['participant'].unique():
        df_sub = df_subset[(df_subset['participant']==part)]

        sample = None

        df_sub = df_sub.reset_index().set_index('indicator')

        df_sub = df_sub[~df_sub.index.isin(['post_tiss', 'post_tiss-1', 'post_tiss-2', 'C1D15'])]

        if 'S' in list(df_sub.index):
            sample = df_sub.loc['S']['analysis_id']
        elif 'C1D1' in list(df_sub.index) and df_sub.iloc[0]['participant']!='MCC095':
            sample = df_sub.loc['C1D1']['analysis_id']
        elif 'tiss' in list(df_sub.index):

            if part == 'OSU02-031':
                sample = 'M031-NR-tiss'
            else:
                sample = df_sub[df_sub['collection_date_dfd']==max(df_sub['collection_date_dfd'])]['analysis_id'].item()

        if sample is not None:
            if sample == 'M095-NR-C1D15':
                df.loc['M095-NR-early_tiss', 'pre-treatment'] = 'yes'
            else:
                df.loc[sample, 'pre-treatment'] = 'yes'

    return df


def truth_table_to_comut_data(df, event_type):

    if 'sample' in list(df.columns):
        df.set_index('sample', inplace=True)

    samples = df.index
    cytobands = df.columns

    df_dict = {}
    df_dict['sample']=[]
    df_dict['category']=[] #cytoband
    for cyt, samp in itertools.product(cytobands, samples):
        if df.loc[samp][cyt] == True:
            df_dict['sample'].append(samp)
            df_dict['category'].append(cyt)

    df_dict['value']=[event_type for q in range(len(df_dict['sample']))]

    return pd.DataFrame(df_dict)


def parse_gtf(gtf_file):

    gtf = pd.read_csv(gtf_file, sep='\t', skiprows=6, header=None)
    gtf = gtf[gtf[2]=='gene']
    gtf[8] = gtf[8].str.extract('((?<=gene_name\W")\S.*?(?="))')
    gtf = gtf[[8, 3, 4, 0]]
    gtf.rename(columns={8:'gene', 3:'start', 4:'end', 0:'chrom'}, inplace=True)
    gtf.set_index('gene', inplace=True)

    return gtf

def list_most_recent_pretreatment_and_posttreat_samples(df):
    '''
    assumes merged annotation file and no index (for comut data)
    '''

    pre_samples = []
    post_samples = []

    # df_subset = df[(df['indicator'].str.contains('tiss'))|(df['indicator'].str.contains('S'))|(df['indicator'].str.contains('C1D1'))]

    for part in df['participant'].unique():
        df_sub = df[(df['participant']==part)]
        df_sub.drop_duplicates('sample', inplace=True)

        pre_sample = None
        post_sample = None

        early_indicators = ['C1D1', 'S', 'tiss-2', 'prim', 'early_tiss', 'met', 'tiss', 'pre_io_prim', 'pre_io_tiss','recur-1', 'recur-2',
                            'tiss-3', 'met-1', 'tiss-1', 'pre_io_met', 'met-2', 'prim-1', 'prim-2']

        df_early = df_sub[df_sub['indicator'].isin(early_indicators)]

        if df_early.shape[0]>0:
            # if df_early.shape[0]>1:
            max_time = max(df_early['collection_date_dfd'])
            # print(df_early[['analysis_id', 'collection_date_dfd', 'participant']])
            if np.isnan(max_time):
                if df_early.shape[0]==1:
                    pre_sample = df_early.set_index('participant').loc[part]['analysis_id']
                else:
                    print('multiple nan samples')
            else:
                pre_sample = df_early[df_early['collection_date_dfd']==max_time]['analysis_id'].item()

        if pre_sample is not None:
            pre_samples.append(pre_sample)


        df_late = df_sub[df_sub['indicator'].isin(['post_tiss', 'EOT', 'C4D1' 'post_tiss-1'])]

        if df_late.shape[0]>0:
            # print(df_late[['analysis_id', 'collection_date_dfd', 'participant']])
            post_sample = df_late[df_late['collection_date_dfd']==max(df_late['collection_date_dfd'])]['analysis_id'].item()

        if post_sample is not None:
            post_samples.append(post_sample)


    return pre_samples, post_samples

def cat_mixcr_clonotypes(S):

    df = pd.concat(
        {
            samp: pd.read_csv(S.loc[samp]['mixcr_all_clones'], sep='\t')
            for samp in S.index
        }
    )

    return df






########
