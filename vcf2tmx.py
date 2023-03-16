"""
A partir de un un archivo vcf y un popmap genera un archivo para análisis de TreeMix
VS
"""

import argparse
import allel
import pandas as pd
import numpy as np
import gzip

#Diccionario de individuos y poblaciones

    popmap = {}
    with open(file_path, 'r') as f:
        for line in f:
            ind, pop = line.strip().split()
            popmap[ind] = pop
    return popmap

#se extraen genotipos

def create_treemix_input(vcf_path, popmap_path, output_path):
    popmap = read_popmap(popmap_path)
    callset = allel.read_vcf(vcf_path, fields=['samples', 'calldata/GT'])
    samples = callset['samples']
    gt = allel.GenotypeArray(callset['calldata/GT'])

    pop_sample_indices = {}
    for i, sample in enumerate(samples):
        if sample in popmap:
            pop = popmap[sample]
            if pop not in pop_sample_indices:
                pop_sample_indices[pop] = []
            pop_sample_indices[pop].append(i)

#poblaciones y alelos 

    pop_gt_counts = {}
    for pop, indices in pop_sample_indices.items():
        pop_gt = gt.take(indices, axis=1)
        pop_gt_counts[pop] = pop_gt.count_alleles()

    with gzip.open(output_path, 'wt') as outfile:
        outfile.write("\t".join(pop_gt_counts.keys()) + "\n")
        for i in range(len(pop_gt_counts[next(iter(pop_gt_counts))])):
            row = []
            for pop, gt_counts in pop_gt_counts.items():
                row.append(f"{gt_counts[i][0]},{gt_counts[i][1]}")
            outfile.write("\t".join(row) + "\n")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert VCF and popmap files to TreeMix input format.')
    parser.add_argument('-v', '--vcf', required=True, help='Input VCF file path.')
    parser.add_argument('-p', '--popmap', required=True, help='Input popmap file path.')
    parser.add_argument('-o', '--output', required=True, help='Output file name for TreeMix input.')

    args = parser.parse_args()

    output_path = f"{args.output}.tmx.gz"
    create_treemix_input(args.vcf, args.popmap, output_path)
