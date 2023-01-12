#!/usr/bin/env python
import argparse
import pandas as pd
import re
import pysam
import pysamstats

def get_coverage(bam, chrom, pos, allele):
    coverage = altcount = refcount = 0;ref=allele[0];alt=allele[1]
    for p in bam.pileup(chrom, pos, pos+1,truncate= True):
        reads=len(p.pileups)
        for pin in p.pileups:
            if not pin.is_del and not pin.is_refskip:
                base=pin.alignment.query_sequence[pin.query_position]
                if base==ref:
                    refcount += 1
                else:
                    if alt.endswith(base):
                        altcount += 1
        coverage=(p.n, reads, refcount, altcount)
    return coverage

def get_rms_mapq(bam, chrom, pos):
    for rec in pysamstats.stat_mapq(bam, chrom=chrom,start=pos,end=pos+1,truncate=True):
            mapq=rec['rms_mapq']
    return mapq

def get_rms_baseq(bam, chrom, pos):
    for rec in pysamstats.stat_baseq(bam, chrom=chrom,start=pos,end=pos+1,truncate=True):
            baseq=rec['rms_baseq']
    return baseq

def cli_interface():
    parser = argparse.ArgumentParser(description='pysam stats.')
    parser.add_argument('--input', dest='input_filename', nargs='+', type=str,
                        help='input batch/sample name', required=True)
    args = parser.parse_args()
    return args

def main():
    args = cli_interface()
    for file in args.input_filename:
        outname=re.sub('$','-coverage.csv',file)
        tumor_dir=re.sub('$','-aPC',file)
        tfilebam=re.sub('$','-ready.bam',tumor_dir)
        tumor_file=("%s/%s" % (tumor_dir,tfilebam))
        normal_file = re.sub('aPC','PB',tumor_file)
        vcf_file = re.sub('$','-ensemble-annotated.vcf.gz',file)
        tumor_bam = pysam.AlignmentFile(tumor_file, "rb")
        normal_bam = pysam.AlignmentFile(normal_file, "rb")
        vcf = pysam.VariantFile(vcf_file)
        results=[]
        for variant in vcf.fetch():
            if variant.id is not None:
                tumor_coverage = get_coverage(tumor_bam, variant.chrom, variant.pos,variant.alleles)
                normal_coverage = get_coverage(normal_bam, variant.chrom, variant.pos,variant.alleles)
                tumor_mapq = get_rms_mapq(tumor_bam, variant.chrom, variant.pos)
                normal_mapq = get_rms_mapq(normal_bam, variant.chrom, variant.pos)
                tumor_baseq = get_rms_baseq(tumor_bam, variant.chrom, variant.pos)
                normal_baseq = get_rms_baseq(normal_bam, variant.chrom, variant.pos)
                tmp_record = [variant.chrom, variant.pos, variant.ref, variant.alleles[1], variant.id,
                              tumor_coverage[0], tumor_coverage[1], tumor_coverage[2], tumor_coverage[3], tumor_mapq,
                              tumor_baseq, normal_coverage[0], normal_coverage[1], normal_coverage[2],
                              normal_coverage[3], normal_mapq, normal_baseq]
                results.append(tmp_record)
        data_frame = pd.DataFrame(results, columns=['chrom', 'pos', 'ref', 'alt', 'rsid', 'tumor_coverage',
                                                    'tumor_PileupRead', 'tumor_rd','tumor_vd', 'tumor_rms_mapq',
                                                    'tumor_rms_baseq', 'normal_coverage', 'normal_PileupRead',
                                                    'normal_rd', 'normal_vd', 'normal_rms_mapq', 'normal_rms_baseq'])

        data_frame.to_csv(outname)
if __name__ == "__main__":
    main()