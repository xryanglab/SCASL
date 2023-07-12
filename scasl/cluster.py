import pandas as pd
import numpy as np
import umap
from sklearn.decomposition import PCA
import argparse
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
from sklearn.cluster import SpectralClustering, KMeans, AgglomerativeClustering
from sklearn.metrics import calinski_harabasz_score
import os

plt.switch_backend('agg')


# sample command
# python cluster.py -m ilc_mat_pca.npy --csv ilc_impute.csv -l ilc_celltype.csv -t facs_gating -o ilc_cluster.csv -n 6 --manifold umap -p 10


def reduce_dim(mat, manifold_alg, pca_dim):
    pca = PCA(n_components=pca_dim)
    embedding = pca.fit_transform(mat).astype(float)
    if manifold_alg == 'umap':
        reducer = umap.UMAP(random_state=42)
    else:
        reducer = TSNE()
    coor = reducer.fit_transform(embedding)
    x = coor[:, 0]
    y = coor[:, 1]
    return embedding, x, y


def make_cluster_labels(embedding, cluster_alg, n_cluster):
    if cluster_alg == 'spectral':
        alg = SpectralClustering(n_clusters=n_cluster, affinity='nearest_neighbors', random_state=7)
    elif cluster_alg == 'agglomerative':
        alg = AgglomerativeClustering(n_clusters=n_cluster)
    else:  # kmeans
        alg = KMeans(n_clusters=n_cluster)
    preds = alg.fit(embedding).labels_
    return preds


def make_label_df(norm_csv, label_file, index_name):
    labels = pd.read_csv(label_file, sep=",")
    labels = labels.set_index(index_name).T
    df_normalized = pd.read_csv(norm_csv).set_index('Site')
    samples = df_normalized.columns.values
    label_samples = labels.columns.values
    diff = set(samples) - set(label_samples)
    labels[list(diff)] = 'undef'
    labels = labels[samples].T
    return labels


def make_sample_df(norm_csv):
    df_normalized = pd.read_csv(norm_csv).set_index('Site')
    samples = df_normalized.columns.values
    labels = pd.DataFrame(samples)
    return labels


def plot_label(labels, label_name, img_path):
    sns.set(style='white')
    plt.figure(figsize=(12, 8))
    assert label_name is not None
    sns.scatterplot(data=labels, x='x', y='y', hue=label_name)
    sns.despine()
    sns.color_palette("Spectral", as_cmap=True)
    reduce_save_path = os.path.join(img_path, 'reduce.png')
    if not os.path.exists(img_path):
        os.makedirs(img_path)
    plt.savefig(reduce_save_path)
    print('dimension reduction plot saved at %s' % (reduce_save_path))
    plt.clf()


def plot_cluster(labels, img_path):
    sns.set(style='white')
    plt.figure(figsize=(10, 8))
    sns.scatterplot(data=labels, x='x', y='y', hue='preds', palette='tab20')
    sns.despine()
    if not os.path.exists(img_path):
        os.makedirs(img_path)
    cluster_save_path = os.path.join(img_path, 'cluster.png')
    plt.savefig(cluster_save_path)
    print('cluster result plot saved at %s' % cluster_save_path)
    plt.clf()


def plot(x, y, labels, preds, label_name, output_path, use_label, img_path):
    print(labels, len(x))
    labels['x'] = x
    labels['y'] = y

    if use_label:
        plot_label(labels, label_name, img_path)

    labels['preds'] = preds
    plot_cluster(labels, img_path)

    labels.to_csv(output_path)


def select(param_range, get_model_pred, metric, features):
    best_score = -1.1
    best_param = None
    best_preds = None
    for param in param_range:
        preds = get_model_pred(param)
        score = metric(features, preds)
        if score > best_score:
            best_score = score
            best_param = param
            best_preds = preds
    return best_param, best_preds, best_score


def cluster(mat_file, norm_csv, label_file, label_name, index_name, output_file, pca_dim, manifold_alg, cluster_alg,
            max_n_cluster, use_label, img_path):
    print('reading file...')
    mat = np.load(mat_file)
    print('done.')

    print('reducing dimensions using PCA and %s...' % (manifold_alg,))
    embedding, x, y = reduce_dim(mat, manifold_alg, pca_dim)
    print('done')

    print('performing cluster with %s algorithm' % (cluster_alg,))

    def make_cluster_pred(n_cluster):
        return make_cluster_labels(embedding, cluster_alg, n_cluster)

    n_cluster, preds, score = select(
        range(2, max_n_cluster + 1),
        make_cluster_pred,
        calinski_harabasz_score, embedding
    )

    print('best n: [%02d] best silhouette score: [%.4f]' % (n_cluster, score))
    preds = make_cluster_labels(embedding, cluster_alg, max_n_cluster)

    num_nan_list = []
    for i in range(max_n_cluster):
        subset = mat[preds == i, :]
        indicator_pos = 2 * np.arange(subset.shape[1] // 2) + 1
        indicator = subset[:, indicator_pos]
        num_nan = np.mean(indicator) * 100
        num_nan_list.append(num_nan)
    print('num nan list: ', num_nan_list)

    print('making label dataframe...')
    if use_label:
        labels = make_label_df(norm_csv, label_file, index_name)
    else:
        labels = make_sample_df(norm_csv)

    print('saving visualization results and saving label file...')
    plot(x, y, labels, preds, label_name, output_file, use_label, img_path)
    print('done.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--mat', default='process_result/normalized_mat.npy', type=str)
    parser.add_argument('--csv', dest='norm_csv', default='process_result/normalized_matrix.csv', type=str)
    parser.add_argument('--no-label', dest='use_label', action='store_false')
    parser.add_argument('-l', '--label', default='data/label.csv', type=str)
    parser.add_argument('-t', '--truth', type=str, default='finaltype')
    parser.add_argument('-o', '--output', type=str, default='process_result/cluster.csv')
    parser.add_argument('-p', '--pca', dest='pca_dim', default=20, type=int)
    parser.add_argument('--manifold', type=str, choices=['umap', 'tsne'], default='umap')
    parser.add_argument('-c', '--cluster', dest='cluster_alg', choices=['spectral', 'agglomerative', 'kmeans'],
                        default='spectral')
    parser.add_argument('-n', '--max_n_cluster', type=int, default=8)
    parser.add_argument('-i', '--index', type=str, default='Run')
    parser.add_argument('--img_path', help='visualization save directory', type=str, default='img')
    args = parser.parse_args()

    cluster(args.mat, args.norm_csv, args.label, args.truth, args.index, args.output, args.pca_dim, args.manifold,
            args.cluster_alg, args.max_n_cluster, args.use_label, args.img_path)
