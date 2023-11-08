# This file aims to turn junction files to junction matrix
# Inputs a series of junction files, outputs a pandas DataFrame

# sample command for unit testing
# python process.py --junc data/junc -o process_result

import numpy as np
import pandas as pd
import os
from .utils import write_json
from tqdm import tqdm
import argparse


def get_sites_from_junction(sites_set, fname):
    df = pd.read_csv(fname, sep='\t', header=None)
    for line in df.values:
        sites_set.add('_'.join([line[0], str(line[1]), str(line[2])]))


def get_sites(junc_dir, output_dir, save=True):
    '''
    return a numpy array of sites extracted from junc files
    '''
    file_list = os.listdir(junc_dir)
    sites_set = set()
    bar = tqdm(file_list)
    bar.set_description('Loading site names')
    for filename in bar:
        if filename.endswith('.junc') or filename.endswith('.tab'):
            get_sites_from_junction(sites_set, os.path.join(junc_dir, filename))
    sites = list(sites_set)
    sites.sort()
    if save:
        write_json(os.path.join(output_dir, 'sites.json'), sites)
    sites_dict = {}
    for i, s in enumerate(sites):
        sites_dict[s] = i
    if save:
        write_json(os.path.join(output_dir, 'sites_dict.json'), sites_dict)
    return sites, sites_dict


def fill_junc_matrix_from_junction(vec, fname, sites_dict):
    df = pd.read_csv(fname, sep='\t', header=None)
    for line in df.values:
        site = '_'.join([line[0], str(line[1]), str(line[2])])
        value = line[4]
        pos = sites_dict[site]
        vec[pos] = value

def fill_junc_matrix_from_junction_star(vec, fname, sites_dict):
    df = pd.read_csv(fname, sep='\t', header=None)
    for line in df.values:
        site = '_'.join([line[0], str(line[1]), str(line[2])])
        value = line[6]
        pos = sites_dict[site]
        vec[pos] = value

def make_junc_matrix(junc_dir, sites, sites_dict, output_dir, save=True):
    file_list = os.listdir(junc_dir)
    num_sites = len(sites)
    num_samples = 0
    mat = []
    samples = []
    bar = tqdm(file_list)
    bar.set_description('Reading and processing junction files')
    for filename in bar:
        if filename.endswith('.junc'):
            num_samples += 1
            sample_name = filename.split('.')[0]
            samples.append(sample_name)
            vec = np.ones(num_sites) * np.nan
            fill_junc_matrix_from_junction(vec, os.path.join(junc_dir, filename), sites_dict)
            mat.append(vec)
        elif filename.endswith('.tab'):
            num_samples += 1
            sample_name = filename.split('.')[0]
            samples.append(sample_name)
            vec = np.ones(num_sites) * np.nan
            fill_junc_matrix_from_junction_star(vec, os.path.join(junc_dir, filename), sites_dict)
            mat.append(vec)
    mat = np.stack(mat)
    if save:
        write_json(os.path.join(output_dir, 'samples.json'), samples)
    np.save(os.path.join(output_dir, 'junc_mat.npy'), mat)
    return mat, samples


def make_dataframe(mat, samples, sites):
    df = pd.DataFrame(mat.T, columns=samples)
    df.insert(0, 'Site', sites)
    return df


def process(junc_dir, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    sites, sites_dict = get_sites(junc_dir, output_dir)
    mat, samples = make_junc_matrix(junc_dir, sites, sites_dict, output_dir)
    return make_dataframe(mat, samples, sites)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--junc', type=str, default='data/junc')
    parser.add_argument('-o', '--output', type=str, default='process_result')
    args = parser.parse_args()
    df = process(args.junc, args.output)
    df.to_csv(os.path.join(args.output, 'junc_matrix.csv'), index=False)
