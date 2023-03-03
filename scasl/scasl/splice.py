from .extract import *
from .process import *
from .filter import *
from .normalize import *
from .cluster import *
import yaml
from easydict import EasyDict
from datetime import datetime
import shutil


def get_config(config_path):
    with open(config_path, 'r') as f:
        cfg = yaml.load(f, yaml.CLoader)
    cfg = EasyDict(cfg)
    cfg = prepare_result_snapshot(cfg)
    shutil.copy(
        os.path.join(config_path),
        os.path.join(cfg.io.output_path)
    )
    return cfg


def prepare_result_snapshot(cfg):
    timestr = datetime.now().strftime("%Y%m%d%H%M%S")
    save_path = os.path.join(cfg.io.output_path, timestr)
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    curdir = os.getcwd()
    filelist = os.listdir(curdir)
    for f in filelist:
        if f.endswith('.py'):
            shutil.copy(
                os.path.join(curdir, f),
                os.path.join(save_path, f)
            )
    cfg.io.output_path = save_path
    return cfg


def run_cluster(cfg):
    # Preprocess
    print('=============Preprocessing=============')
    extract(cfg.process.bam, cfg.process.junction, cfg.process.lc)
    df = process(cfg.process.junction, cfg.io.output_path)
    junc_path = os.path.join(cfg.io.output_path, 'junc_matrix.csv')
    df.to_csv(junc_path, index=False)

    # Filter
    print('=============Filtering=============')
    img_path = os.path.join(cfg.io.output_path, 'img')
    if not os.path.exists(img_path):
        os.makedirs(img_path)
    df_start, df_end = filter(
        junc_path,
        cfg.threshold.sites_initial, cfg.threshold.runs_initial,
        cfg.threshold.sites_quality, cfg.threshold.runs_quality,
        cfg.io.log_scale_histogram, img_path
    )
    print('saving...')
    filter_path = os.path.join(cfg.io.output_path, 'filtered_matrix')
    df_start.to_csv(filter_path + '_start.csv', index=False)
    df_end.to_csv(filter_path + '_end.csv', index=False)
    print('done.')

    # Normalize and Impute
    print('=============Normalization & Imputation=============')
    df_final, mat = normalize(filter_path, cfg.impute.num_iteration, cfg.impute.knn)

    print('saving files...')
    normalize_path = os.path.join(cfg.io.output_path, 'normalized_matrix.csv')
    mat_path = os.path.join(cfg.io.output_path, 'normalized_matrix.npy')
    df_final.to_csv(normalize_path)
    np.save(mat_path, mat)
    print('done.')

    # Cluster
    print('=============Cluster & Visualize=============')
    result_path = os.path.join(cfg.io.output_path, 'cluster_result.csv')
    cluster(
        mat_path, normalize_path,
        cfg.io.label_file, cfg.io.truth_column, cfg.io.run_column,
        result_path,
        cfg.cluster.pca_dimension,
        cfg.cluster.dimension_reduction_method,
        cfg.cluster.cluster_algorithm,
        cfg.cluster.num_cluster,
        cfg.io.use_label, img_path
    )


class SCASL:
    def __init__(self, cfg):
        self.cfg = cfg

    def fit(self):
        run_cluster(self.cfg)