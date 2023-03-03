# This file is used to extract junction reads to make junction matrix

import os
import argparse


def extract(bam_dir, junc_dir, lc_path):
    if not os.path.exists(junc_dir):
        os.makedirs(junc_dir)
    os.system('ls %s/*.bam > bam_path.txt' % (bam_dir, ))
    shell_cmd = 'cat bam_path.txt |while read id\n'\
    'do\n'\
    'file=$(basename $id )\n'\
    'sample=${file%%%%.*}\n'\
    '\techo Converting $id to $sample.junc\n'\
    '\tsh %s/scripts/bam2junc.sh  $id $sample.junc\n'\
    'done' % (lc_path,)
    'ls *.junc > all_juncfiles.txt'
    print(shell_cmd)
    os.system(shell_cmd)
    os.system('mv *.junc %s' % (junc_dir, ))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-j', '--junc', type=str, default='junction', help='junction file directory')
    parser.add_argument('--bam', type=str, required=True, help='bam file directory')
    parser.add_argument('--lc', type=str, required=True, help='leafcutter path')
    args = parser.parse_args()
    extract(args.bam, args.junc, args.lc)
