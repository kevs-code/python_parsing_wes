#!/usr/bin/env python
import os
import sys
import pysam
import pysamstats
#import matplotlib.pyplot as plt
import re

def get_coverage(bam, chrom, pos, allele):
    coverage = altcount = refcount = altcount2 = count2 = 0
    ref=allele[0];alt=allele[1]
    for p in bam.pileup(chrom, pos, pos+1,truncate= True):
        print ("\ncoverage at base %s = %s" %
           (p.pos, p.n))
        reads=len(p.pileups)
        #print(p)
        for pin in p.pileups:
            #print (pin)
            count2 +=1
            if not pin.is_del and not pin.is_refskip:
                base=pin.alignment.query_sequence[pin.query_position]
                if base==ref:
                    print ('\tbase %s ref %s alt %s' % (base,ref,alt))
                    refcount += 1
                else:
                    altcount += 1
                if alt.endswith(base):
                    altcount2 += 1
        coverage = p.n
        print(coverage, refcount, altcount, altcount2 ,count2)
        coverage=(p.n, reads, refcount, altcount)
    return coverage

def get_coverage_ext(bam, chrom, pos):
    for rec in pysamstats.stat_coverage_ext(bam, chrom=chrom,start=pos,end=pos+1,truncate=True):
            print (rec)

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

os.chdir("/home/kevin/mrd_wes/results/2018-05-09/samples")
tumor_file = "MRD-16-028-aPC/MRD-16-028-aPC-ready.bam"
normal_file = "MRD-16-028-PB/MRD-16-028-PB-ready.bam"
vcf_file = "../MRD-16-028-ensemble-annotated.vcf.gz"

tumor_bam = pysam.AlignmentFile(tumor_file, "rb")
normal_bam = pysam.AlignmentFile(normal_file, "rb")
vcf = pysam.VariantFile(vcf_file)
results=[]
for variant in vcf.fetch():
    if variant.id is not None:
        tumor_coverage = get_coverage(tumor_bam, variant.chrom, variant.pos,variant.alleles)
        #get_coverage_ext(tumor_bam, variant.chrom, variant.pos)
        normal_coverage = get_coverage(normal_bam, variant.chrom, variant.pos)
        tumor_mapq = get_rms_mapq(tumor_bam, variant.chrom, variant.pos)
        normal_mapq = get_rms_mapq(normal_bam, variant.chrom, variant.pos)
        tumor_baseq = get_rms_baseq(tumor_bam, variant.chrom, variant.pos)
        normal_baseq = get_rms_baseq(normal_bam, variant.chrom, variant.pos)
        tmp_record = [variant.chrom, variant.pos, variant.ref, variant.alleles[1], variant.id, tumor_coverage,
                      tumor_mapq, tumor_baseq, normal_coverage, normal_mapq, normal_baseq]
        results.append(tmp_record)