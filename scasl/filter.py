import numpy as np
import pandas as pd
import argparse
import seaborn as sns
import matplotlib.pyplot as plt
import os
from tqdm import tqdm

plt.switch_backend('agg')


def split_start_end(name):
    splits = name.split('_')
    start = '_'.join(splits[:2])
    end = '_'.join([splits[0], splits[2]])
    return start, end


def make_start_end_df(df):
    """
    find start & end from the site name
    and add two columns in the dataframe for them
    """
    sites = df['Site']
    starts_and_ends = list(map(split_start_end, sites))
    df['start'] = list(map(lambda x: x[0], starts_and_ends))
    df['end'] = list(map(lambda x: x[1], starts_and_ends))
    return df


def repeat_filter(df, by):
    """
    return a copy of the original dataframe
    removing duplicated site starts / ends
    """
    assert by in ['start', 'end']
    start_end_df = make_start_end_df(df)
    filtered_df = start_end_df[start_end_df.duplicated(by, keep=False)]
    return filtered_df


def draw_hist(count, split='site', log_scale=True, img_path='img'):
    sns.set()
    plt.figure(figsize=(10, 6))
    sns.histplot(np.log(count) if log_scale else count, bins=30, kde=True)
    if not os.path.exists('img'):
        os.makedirs('img')
    plt.xlabel('log value of the number of non-NaN data items' if log_scale else 'number of non-NaN data items')
    plt.savefig(os.path.join(img_path, '%s_hist.png' % (split,)))
    print('the %s histogram is saved at %s/%s_hist.png' % (split, img_path, split))
    print('the descriptions of the non-NaN data of sites are shown below')
    print(count.describe(percentiles=np.arange(0, 100, 10) / 100.))
    plt.clf()


def qc_filter_sites(df_repeat, thres, log_scale=True, img_path='img'):
    site_count = df_repeat.count(axis=1)
    draw_hist(site_count, 'site', log_scale, img_path)
    df_site = df_repeat.copy()
    df_site['site_count'] = site_count.values
    df_site_filtered = df_site[df_site['site_count'] > thres]
    return df_site_filtered


def qc_filter_samples(df_site_filtered, thres, log_scale=True, img_path='img'):
    sample_count = df_site_filtered.count(axis=0)
    draw_hist(sample_count, 'sample', log_scale, img_path)
    df_sample = df_site_filtered.copy()
    df_sample.loc['sample_count'] = sample_count.values
    df_sample_filtered = df_sample.loc[:, df_sample.loc['sample_count'] > thres]
    return df_sample_filtered


def thres_filter(df, samples_ps, sites_ps, by='start'):
    groups = df.groupby(by=by)
    result = []
    for s in tqdm(groups):
        r = s[1].copy().iloc[:, 1:-2]
        sums = r.sum(axis=0)
        counts = r.count(axis=1)
        r[counts < samples_ps] = np.nan
        r.loc[:, sums < sites_ps] = np.nan
        result.append(r)
    result = pd.concat(result)
    result.insert(0, 'Site', df['Site'])
    result[['start', 'end']] = df[['start', 'end']]
    return result


def filter(junc_mat, samples_ps, sites_ps, sites_thres, samples_thres, log_scale, img_path):
    print('reading file...')
    df = pd.read_csv(junc_mat, index_col=False)
    # import pdb; pdb.set_trace()
    print('done.')
    print('executing repeat and initial threshold filter...')
    start_df = repeat_filter(df, 'start')
    end_df = repeat_filter(df, 'end')
    start_df = thres_filter(start_df, samples_ps, sites_ps, by='start')
    end_df = thres_filter(end_df, samples_ps, sites_ps, by='end')
    print('done')
    print('executing sites quality filter by threshold')
    start_df = qc_filter_sites(start_df, sites_thres, log_scale=log_scale, img_path=img_path)
    end_df = qc_filter_sites(end_df, sites_thres, log_scale=log_scale, img_path=img_path)
    print('done.')

    print('remove the duplicated site starts and ends...')
    start_df = repeat_filter(start_df.iloc[:-1, :-3], 'start')
    end_df = repeat_filter(end_df.iloc[:-1, :-3], 'end')
    print('done.')

    print('executing sample quality filter...')
    # sample filter needs the combined dataframe
    df_repeat = pd.concat([start_df, end_df])
    df_sample_filtered = qc_filter_samples(df_repeat, samples_thres, log_scale=log_scale, img_path=img_path)
    print('done.')

    # use the filtered sites to selec from the original dataframe
    start_df = start_df[df_sample_filtered.columns]
    end_df = end_df[df_sample_filtered.columns]

    # executing the final repeat filter
    start_df = repeat_filter(start_df, 'start').sort_values(by='start')
    end_df = repeat_filter(end_df, 'end').sort_values(by='end')

    return start_df, end_df


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--junc', default='process_result/junc_matrix.csv', type=str)
    parser.add_argument('-o', '--output', default='process_result/filtered_matrix', type=str)
    parser.add_argument('--st', dest='sites_ps', default=1, type=int)
    parser.add_argument('--sp', dest='samples_ps', default=1, type=int)
    parser.add_argument('--sites', dest='sites_thres', default=3, type=int)
    parser.add_argument('--samples', dest='samples_thres', default=5, type=int)
    parser.add_argument('--log', action='store_true',
                        help='whether to perform logorithm scale on the counts for the histogram plotting')
    parser.add_argument('--img_path', help='visualization save directory', type=str, default='img')
    args = parser.parse_args()
    dfs, dfe = filter(args.junc, args.samples_ps, args.sites_ps, args.sites_thres, args.samples_thres, args.log, args.img_path)
    print('saving...')
    dfs.to_csv(args.output + '_start.csv', index=False)
    dfe.to_csv(args.output + '_end.csv', index=False)
    print('done.')
