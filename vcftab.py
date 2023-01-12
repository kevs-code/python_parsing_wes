#!/usr/bin/env python
import argparse
from collections import defaultdict
import numpy as np
import pandas as pd
import vcf

def create_mutation_index(filenames):
    varindex = defaultdict(list)

    for vf in filenames:
        vcf_reader = vcf.Reader(filename=vf)
        for record in vcf_reader:
            varindex[(record.CHROM, record.POS, record.REF,str(record.ALT))].append(record)

    return varindex


def get_samples_by_record(records):
    samples = set()
    for r in records:
        for rs in r.samples:
            if rs.sample.endswith('aPC'):
                samples.add(rs.sample)

    return samples

def get_genes_by_record(records):
    genes = set()
    for r in records:
        for ann in r.INFO['ANN']:
            ann_fields = ann.split('|')[3]
            genes.add(ann_fields)
    return genes

def get_mean_dp_by_record(records):
    DP = []
    for r in records:
        for rs in r.samples:
            if rs.sample.endswith('aPC'):
                DP.append(r.genotype(rs.sample)['DP'])

    return np.mean(DP)

def get_mean_sample_field_by_record(records, field):
    val = []
    for r in records:
        for rs in r.samples:
            if rs.sample.endswith('aPC'):
                try:
                    val.append(rs[field])
                except AttributeError:
                    None

    if len(val) > 0:
        return np.mean(val)
    else:
        return 'NA'


def cli_interface():
    parser = argparse.ArgumentParser(description='VCF to Table.')
    parser.add_argument('--input', dest='input_filename', nargs='+', type=str, help='input VCF files', required=True)
    parser.add_argument('--output', dest='output_filename', type=str, help='output CSV files', required=True)
    args = parser.parse_args()
    return args


def main():
    args = cli_interface()

    # creating the mutation index
    varindex = create_mutation_index(args.input_filename)

    results = []
    # processing the index
    for mut, records in varindex.items():
        # building sample info
        samples = get_samples_by_record(records)
        genes = get_genes_by_record(records)
        dp_info = get_mean_dp_by_record(records)
        mq_info = get_mean_sample_field_by_record(records, 'MQ')
        mmq_info = get_mean_sample_field_by_record(records, 'MMQ')
        for g in genes:
            tmp_record = [mut[0], mut[1], mut[2], mut[3], dp_info, mq_info, mmq_info, len(samples), ";".join(samples), g]
            results.append(tmp_record)

    # converting to dataframe and then to csv
    data_frame = pd.DataFrame(results, columns=['chrom', 'pos', 'ref','alt', 'dp', 'mq', 'mmq', 'num_samples', 'samples', 'Gene_Name' ])

    data_frame.to_csv(args.output_filename)

if __name__ == "__main__":
    main()
