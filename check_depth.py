#!/usr/bin/env python

import sys
import pysam
import pysamstats

def get_coverage(bam, chrom, pos):
    coverage = 0
    for p in bam.pileup(chrom, pos, pos+1,truncate= True):
        coverage = p.n

    return coverage

def main():
    tumor_file = sys.argv[1]
    normal_file = sys.argv[2]
    vcf_file = sys.argv[3]

    tumor_bam = pysam.AlignmentFile(tumor_file, "rb")
    normal_bam = pysam.AlignmentFile(normal_file, "rb")
    vcf = pysam.VariantFile(vcf_file)

    for variant in vcf.fetch():
        if variant.id is not None:
            tumor_coverage = get_coverage(tumor_bam, variant.chrom, variant.pos)
            normal_coverage = get_coverage(normal_bam, variant.chrom, variant.pos)

            print(variant.chrom, variant.pos, variant.id, tumor_coverage, normal_coverage)


if __name__ == "__main__":
    main()
