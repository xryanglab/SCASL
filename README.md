# SCASL: single-cell clustering based on alternative splicing landscapes #
SCASL is a strategy of cell clustering by systematically assessing the landscapes of single-cell RNA splicing. SCASL is mainly used to 1. identify novel single cell cluster based on the AS level, 2. reveal the transition relationship between clusters and the differential splicing pattern, 3. identify important differential splicing events.

If you want to use this method or our result in your research, please cite our paper, the detailed introduction and application of the algorithm can also be found in the paper.

## Version
1.0.0 (We will continue to supplement the functions and application range of scasl)

## Author
Xianke Xiang, Xuerui Yang.

## Getting start ##

### Dependencies and requirements
The packages that SCASL depends on and the versions corresponding to the packages are organized in `requirements.txt`. The package has been tested on conda 4.5.11 and is platform independent (tested on Windows, macOS and Linux). 

### Environment Setup
```bash
> conda create -n scasl -c conda-forge scikit-learn python=3.9
> conda activate scasl
> conda install -c conda-forge pandas pyyaml seaborn tqdm easydict umap-learn
> conda install -c davidaknowles r-leafcutter
```

It should be noted that if leafcutter cannot be downloaded successfully with conda, you can also consider using source code.
```bash
> git clone https://github.com/davidaknowles/leafcutter
```

## Run
```bash
> conda activate scasl
> python main.py -y <your_config_file_path>
```

### Parameter settings
By default, this method uses the bam file (any mapping method can be adopted) as the initial input and can directly output the clustering results. At the same time, it also supports users to use intermediate files as input to run the py files in the scasl folder step by step to get the results they want. 

The following parameters can be adjusted directly in `srr.yaml` to use scasl:

- **bam**: bam file directory.
- **lc**: directory of leafcutter.
- **sites_initial**: quality control of the initial junction matrix, minimum number of expressing cells per site.
- **runs_initial**: quality control of the initial junction matrix, minimum number of expressing sites per cell.
- **sites_quality**: quality control of the AS matrix, minimum number of expressing cells per AS group.
- **runs_quality**: quality control of the AS matrix, minimum number of expressing sites (sum of AS group) per cell.
- **log_scale_histogram**: whether to output the quality distribution histogram to help determine the quality control standard, logical parameter (T or F).
- **use_label**: whether to use known labels when clustering, logical parameter (T or F).
- **label_file**: directory of label file.
- **truth_column**: the column name of the label to use in the label file.
- **run_column**: the column names of the cells (corresponding to the bam file) in the label file.
- **output_path**: the directory of all output files.
- **num_iteration**: the number of iterations during weighted KNN imputation.
- **knn**: the number of adjacent cells used in weighted KNN imputation.
- **pca_dimension**: dimension of PCA when clustering.
- **dimension_reduction_method**: methods for clustering visualization ('umap', 'tsne').
- **cluster_algorithm**: 'spectral' by default. ('spectral', 'agglomerative', 'kmeans'. )
- **num_cluster**: number of clusters.

There are also some other parameters that users can adjust by themselves to optimize the clustering effect according to their own data, such as setting **max_n_cluster** in `cluster.py` to output clustering scores to assist in judging the most appropriate number of clusters.

## Result
The final output file is a cluster label file, **preds** represents the predicted AS cluster label. At the same time, many intermediate files will be generated (such as AS probability matrix, NA position information), and users can extract intermediate files as needed.

As an example of configuration file, `configs/srr.yaml` and [bam](https://drive.google.com/drive/folders/1sFBoileBgYH46QiW6mohR82fr4DUhzGJ?usp=sharing) provides a minimized version of bam data for scRNA-seq. You can also choose to use more cells (from TNBC-2) in `data/junction` for testing.

## License
SCASL is licensed under the Apache License 2.0.
