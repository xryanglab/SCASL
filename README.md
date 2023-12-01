# SCASL: single-cell clustering based on alternative splicing landscapes #
SCASL is a strategy of cell clustering by systematically assessing the landscapes of single-cell RNA splicing. SCASL is mainly used to: 1) Extraction of AS information from scRNA-seq data and imputation of missing values, 2) identify novel single cell cluster based on the AS level, 3) reveal the transition relationship between clusters and the differential splicing pattern, 4) identify important differential splicing events.

If you want to use this method or our result in your research, please cite our paper, the detailed introduction and application of the algorithm can also be found in the paper.

## Version
1.0.0 (We will continue to supplement the functions and application range of scasl)

## Author
Xianke Xiang, Xuerui Yang.

## Getting start ##

### Dependencies and requirements
The packages that SCASL depends on and the versions corresponding to the packages are organized in `requirements.txt`. The package has been tested on conda 4.10.3 and is platform independent (tested on Windows, macOS and Linux). 

### Environment Setup
```bash
> conda create -n scasl -c conda-forge scikit-learn python=3.9
> source activate scasl
> conda install -c conda-forge pandas pyyaml seaborn tqdm easydict umap-learn
> conda install -c davidaknowles r-leafcutter
> conda install -c bioconda samtools -y
```

**Tips** 

1. If you encounter difficulties downloading Leafcutter using conda, you have the option to use the source code instead and add the path to the source code to the environment variable.
```bash
> git clone https://github.com/davidaknowles/leafcutter
> export PATH="YOUR_PATH_OF_LEAFCUTTER:$PATH"
```

2. In addition to using Leafcutter, users also have the option to utilize the 'SJ.out.tab' files generated automatically from the STAR mapping pipeline as junction files. To modify the **junction path**, simply navigate to the location of the "SJ.out.tab" file in the `configs/srr.yaml` file.

It usually takes less than 15 minutes to complete the environment configuration. 

## Run
```bash
> source activate scasl
> python main.py -y configs/srr.yaml
```

### Parameter settings
By default, this method uses the **bam** file (any mapping method can be adopted) or **junction** file as the **initial input** and can directly output the clustering results. At the same time, it also supports users to use intermediate files as input to run the py files in the scasl folder step by step to get the results they want. 

The following parameters can be adjusted directly in `configs/srr.yaml` to use scasl:

- **bam**: bam file directory.
- **lc**: directory of leafcutter.
- **junction**: junction file directory. If no junction directory is provided, a junction directory will be automatically generated in the current location.
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
The final **output** file is a cluster label file, **preds** represents the predicted **AS cluster** label. At the same time, many intermediate files will be generated (such as AS probability matrix, NA position information), and users can extract intermediate files as needed.

Elucidation of intermediate files:
- **site_hist.png**: This figure is used to assist in the selection of filter parameters. The horizontal axis denotes the count of non-NA samples at an individual site, while the vertical axis represents the number of such sites.
- **site_hist.png**: This figure is used to assist in the selection of filter parameters. The horizontal axis denotes the count of non-NA sites at an individual sample, while the vertical axis represents the number of such samples.
- **reduce.png**: This figure represents the outcome of the final dimensionality reduction visualization, where the coloring is based on the ground truth provided by the user.
- **cluster.png**: This figure represents the outcome of the final dimensionality reduction visualization, where the coloring is based on the clustering results of SCASL.
- **junc_matrix.csv,junc_mat.npy**: Junction reads matrix integrated from the single-cell file. In this matrix, the rows correspond to sites, while the columns represent cells.
- **filtered_matrix_start,filtered_matrix_end**: Junction reads matrix after quality filtering and AS grouping, the rows correspond to sites, while the columns represent cells. Start represents grouping by upstream, end represents grouping by downstream.
- **normalized_matrix**: The final AS probability matrix, the rows correspond to sites, while the columns represent cells.

As an example of configuration file, `configs/srr.yaml` and [bam](https://drive.google.com/drive/folders/1sFBoileBgYH46QiW6mohR82fr4DUhzGJ?usp=sharing) provides a minimized version of bam data for scRNA-seq. You can also choose to use the demo of the intermediate files (from TNBC-2) in `data/junction` for testing.

The total test time of the demo files is expected to be less than 3 minutes, and the results are stored in `process_result_demo`. However, github is not suitable for uploading large amounts of data, so the test files used as demo are very few and are only used to show the operation and speed of the software. Due to the randomness of the interpolation, dimensionality reduction and unsupervised clustering processes, there may be some differences between the running results, but the results are generally consistent when the parameters are completely consistent.

## License
SCASL is licensed under the Apache License 2.0.
