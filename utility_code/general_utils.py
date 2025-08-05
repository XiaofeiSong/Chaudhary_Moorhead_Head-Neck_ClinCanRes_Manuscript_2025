import numpy as np
import pandas as pd
import os
import itertools



def list_most_recent_samples(df, indicators, screening_overwrite=False, tissue_priority=False):

    samples = []

    if df.index.name=='analysis_id':
        df['analysis_id'] = df.index

    for part in df['participant'].unique():
        df_sub = df[(df['participant']==part) & (df['indicator'].isin(indicators))]
        df_sub = df_sub.drop_duplicates('analysis_id')

        if tissue_priority and 'tissue' in df_sub['pdb_original_material_type']:
            df_sub = df_sub[df_sub['pdb_original_material_type']=='tissue']

        sample = None

        if df_sub.shape[0]>0:
            max_time = max(df_sub['collection_date_dfd'])
            if np.isnan(max_time):
                if df_sub.shape[0]==1:
                    sample = df_sub.set_index('participant').loc[part]['analysis_id']
                else:
                    print('multiple nan time samples')
            else:
                if df_sub[df_sub['collection_date_dfd']==max_time]['analysis_id'].shape[0]>1:
                    sample = df_sub[(df_sub['indicator']=='prim')|(df_sub['indicator']=='tiss-1')]['analysis_id'].item()
                else:
                    sample = df_sub[df_sub['collection_date_dfd']==max_time]['analysis_id'].item()

        if sample is not None:
            if screening_overwrite:
                if 'S' in df_sub['indicator'].to_list():
                    sample = df_sub[df_sub['indicator']=='S']['analysis_id'].item()

            samples.append(sample)


    return samples



def list_most_recent_pretreatment_samples(df, screening_overwrite=False, tissue_priority=False):

    indicators = ['C1D1', 'S', 'tiss-2', 'prim', 'early_tiss', 'met', 'tiss', 'pre_io_prim', 'pre_io_tiss','recur-1', 'recur-2',
                'tiss-3', 'met-1', 'tiss-1', 'pre_io_met', 'met-2', 'prim-1', 'prim-2']

    return list_most_recent_samples(df, indicators, screening_overwrite=screening_overwrite, tissue_priority=tissue_priority)


def list_most_recent_tissue_samples(df):

    indicators = ['tiss-2', 'prim', 'early_tiss', 'met', 'tiss', 'pre_io_prim', 'pre_io_tiss','recur-1', 'recur-2',
                'tiss-3', 'met-1', 'tiss-1', 'pre_io_met', 'met-2', 'prim-1', 'prim-2']

    return list_most_recent_samples(df, indicators)


def list_most_recent_posttreatment_samples(df):

    indicators = ['post_tiss', 'EOT', 'C4D1', 'post_tiss-1']

    return list_most_recent_samples(df, indicators)


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







########
