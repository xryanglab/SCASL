import numpy as np
import pandas as pd
from tqdm import trange
import argparse


def to_prob_slice(s):
    sums = s.sum(axis=0)
    r = s / sums
    return r


def get_begin_pos(starts):
    eq = list(map(lambda i: starts[i] != starts[i - 1] if i > 0 else True, range(len(starts))))
    begin_pos = np.where(eq)[0].astype(int)
    return begin_pos


def to_prob(df, groupby):
    sums = df.groupby(groupby).sum(min_count=1)
    sums = pd.merge(df, sums, how='left', on=groupby)
    sums = sums.drop(columns=['start', 'end'])
    df_copy = df.copy()
    df_copy = df_copy.drop(columns=['start', 'end'])
    num_samples = df_copy.shape[1]
    sums = sums.iloc[:, num_samples:]
    probs = df_copy.values / sums.values
    mask = ((np.isnan(probs)) & (~np.isnan(sums)))
    probs[mask] = 0
    df_result = df.copy()
    df_result.iloc[:, :-2] = probs
    df_result = df_result.sort_values(by=groupby)
    return df_result


def norm_only(df_path, groupby):
    print(f'reading data from {df_path}...')
    df = pd.read_csv(df_path + '_' + groupby + '.csv', index_col=False)
    df = df.set_index('Site')
    df_prob = to_prob(df, groupby=groupby)
    df_prob = df_prob.drop(['start', 'end'], axis=1)
    print('done.')
    return df_prob


def fill_na(df_prob):
    df_fillna = df_prob.copy()
    for i in trange(len(df_prob)):
        mean = np.nanmean(df_prob.iloc[i, :].values)
        if np.isnan(mean):
            mean = 0.5
        df_fillna.iloc[i, :] = df_fillna.iloc[i, :].fillna(mean)
        na_map = df_prob.iloc[i, :].isna().astype(int)
        random_noise = na_map * np.random.randn(len(df_prob.iloc[i, :])) * 0.01
        df_fillna.iloc[i, :] += random_noise
    df_fillna[df_fillna > 1.] = 1.
    df_fillna[df_fillna < 0] = 0
    return df_fillna


def knn_impute(df_fillna, df_prob, k=5):
    array_fillna = df_fillna.values.T
    m, n = array_fillna.shape
    d = np.zeros((m, m))
    for i in trange(m - 1):
        for j in range(i, m):
            d[i][j] = d[j][i] = np.sqrt(np.sum((array_fillna[i] - array_fillna[j]) ** 2))

    knn = np.zeros((m, k))
    knn_d = np.zeros((m, k))
    for i in range(m):
        knn[i] = np.argsort(d[i])[1:k+1]
        knn_d[i] = np.sort(d[i])[1:k+1]
    knn = knn.astype(int)
    knn_d_inv = np.exp(-2 * knn_d / knn_d.std())
    knn_weights = knn_d_inv / (np.sum(knn_d_inv, axis=1).reshape(m, 1).repeat(k, axis=1))

    knn_val = np.zeros((m, n))
    for i in range(m):
        knn_val[i] = np.average(array_fillna[knn[i]], axis=0, weights=knn_weights[i])
    
    df_knn = df_prob.copy()
    for i in trange(m): 
        na_map = df_prob.iloc[:, i].isna().astype(int)
        df_knn.iloc[:, i] = df_knn.iloc[:, i].fillna(0)
        df_knn.iloc[:, i] += knn_val[i] * na_map.values

    return df_knn


def normalize(df_path, num_iter, k):
    dfs = norm_only(df_path, 'start')
    dfe = norm_only(df_path, 'end')
    df_prob = pd.concat([dfs, dfe])

    print('encoding matrix...')

    df_fillna = fill_na(df_prob)
    array_fillna = df_fillna.values.T
    m, n = array_fillna.shape

    for i in range(num_iter):
        print(f'\n=================== Iter {i + 1} ===================')
        df_fillna = knn_impute(df_fillna, df_prob, k)

    indicator = df_prob.isna().astype(int)
    mat = np.zeros((m, n, 2))
    mat[:, :, 0] = df_fillna.values.T * 1
    mat[:, :, 1] = indicator.values.T * 1
    mat = mat.reshape((m, n * 2))

    print('done.')

    return df_fillna, mat


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--filter', default='process_result/filtered_matrix', type=str)
    parser.add_argument('-o', '--output', default='process_result/normalized_matrix.csv', type=str)
    parser.add_argument('-m', '--mat', default='process_result/normalized_mat.npy', type=str)
    parser.add_argument('-k', default=5, type=int)
    parser.add_argument('-n', '--num-iter', dest='num_iter', default=3, type=int)
    args = parser.parse_args()

    df_final, mat = normalize(args.filter, args.num_iter, args.k)
    print('saving files...')
    df_final.to_csv(args.output)
    np.save(args.mat, mat)
    print('done.')
