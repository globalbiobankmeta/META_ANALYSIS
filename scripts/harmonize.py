#!/usr/bin/env python3
import datetime
import argparse
import json
import gzip
from collections import namedtuple, defaultdict
import sys
import math
from scipy.stats import chi2
import scipy.stats
import numpy
from typing import Dict, Tuple, List
import subprocess
from collections import deque
import re

flip = {"A":"T","C":"G","T":"A","G":"C"}

def check_eff_field(field):
    if field.lower() in ["beta","or"]:
        return field.lower()
    else:
        raise Exception("effect_type must be beta or OR")

def flip_strand( allele):
    return "".join([ flip[a] for a in allele])

def is_symmetric(a1, a2):
    return (a1=="A" and a2=="T") or (a1=="T" and a2=="A") or (a1=="C" and a2=="G") or (a1=="G" and a2=="C")

class Variant():

    def __init__(self, chr, pos, ref, alt, af):
        self.chr = chr
        self.pos = int(float(pos))
        self.ref = ref.strip().upper()
        self.alt = alt.strip().upper()
        self.af = af

    def __eq__(self, other):

        return self.chr == other.chr and self.pos == other.pos and self.ref == other.ref and self.alt == other.alt

    def __lt__(self, other):

        return (  (self.chr==other.chr and self.pos<other.pos)
                  or (self.chr < other.chr)
               )

    def is_equal(self, other:'VariantData') -> bool:
        """
            Checks if this VariantData is the same variant (possibly different strand or ordering of alleles)
            returns: true if the same false if not
        """

        if (self.chr == other.chr and self.pos == other.pos):
            flip_ref =  flip_strand(other.ref)
            flip_alt =  flip_strand(other.alt)

            if self.ref== other.ref and self.alt == other.alt :
                return True

            if is_symmetric( other.ref, other.alt ):
                ## never strandflip symmetrics. Assumed to be aligned.
                if self.ref == other.ref and self.alt == other.alt:
                    return True
                elif self.ref == other.alt and self.alt == other.ref:
                    return True

            elif (self.ref == other.alt and self.alt == other.ref) :
                return True
            elif (self.ref == flip_ref and self.alt==flip_alt):
                return True
            elif (self.ref == flip_alt and self.alt==flip_ref):
                return True

        return False
    

class VariantData(Variant):

    def __init__(self, chr, pos, ref, alt, af, beta, se, pval, extra_cols=[]):
        self.chr = chr
        self.pos = int(float(pos))
        self.ref = ref.strip().upper()
        self.alt = alt.strip().upper()
        self.af = float(af)
        self.beta = float(beta)
        self.se = float(se)
        self.pval = float(pval)
        self.extra_cols = extra_cols

    def equalize_to(self, other:'Variant') -> bool:
        """
            Checks if this VariantData is the same variant as given other variant (possibly different strand or ordering of alleles)
            If it is, changes this variant's alleles, beta and af accordingly
            returns: true if the same (flips effect direction, ref/alt alleles and af if necessary) or false if not the same variant
        """

        if (self.chr == other.chr and self.pos == other.pos):
            flip_ref =  flip_strand(other.ref)
            flip_alt =  flip_strand(other.alt)

            if self.ref== other.ref and self.alt == other.alt :
                return True

            if is_symmetric( other.ref, other.alt ):
                ## never strandflip symmetrics. Assumed to be aligned.
                if self.ref == other.ref and self.alt == other.alt:
                    return True
                elif self.ref == other.alt and self.alt == other.ref:
                    self.beta = -1 * self.beta if self.beta is not None else None
                    self.af = 1 - self.af if self.af is not None else None
                    t = self.alt
                    self.alt = self.ref
                    self.ref = t
                    return True

            elif (self.ref == other.alt and self.alt == other.ref) :
                self.beta = -1 * self.beta if self.beta is not None else None
                self.af = 1 - self.af if self.af is not None else None
                t = self.alt
                self.alt = self.ref
                self.ref = t
                return True
            elif (self.ref == flip_ref and self.alt==flip_alt):
                self.ref = flip_strand(self.ref)
                self.alt = flip_strand(self.alt)
                return True
            elif (self.ref == flip_alt and self.alt==flip_ref):
                self.beta = -1 * self.beta if self.beta is not None else None
                self.af = 1 - self.af if self.af is not None else None
                self.ref =flip_strand(self.alt)
                self.alt = flip_strand(self.ref)
                return True

        return False

def harmonize(file_in, file_ref):

    fp_ref = gzip.open(file_ref, 'rt')
    ref_has_lines = True
    ref_chr = 1
    ref_pos = 0
    ref_h_idx = {h:i for i,h in enumerate(fp_ref.readline().strip().split('\t'))}
    
    with gzip.open(file_in, 'rt') as f:
        h_idx = {h:i for i,h in enumerate(f.readline().strip().split('\t'))}
        print('#CHR\tPOS\tAllele1\tAllele2\tAF_Allele2\tBETA\tSE\tp.value\tAF_gnomAD')
        for line in f:
            s = line.strip().split('\t')
            var = VariantData(s[h_idx['CHR']].replace('chr', '').replace('X', '23'), s[h_idx['POS']],
                              s[h_idx['Allele1']], s[h_idx['Allele2']],
                              s[h_idx['AF_Allele2']], s[h_idx['BETA']],
                              s[h_idx['SE']], s[h_idx['p.value']])
            ref_vars = []
            while ref_has_lines and int(ref_chr) < int(var.chr) or (int(ref_chr) == int(var.chr) and ref_pos < var.pos):
                ref_line = fp_ref.readline().strip().split('\t')
                try:
                    ref_chr = ref_line[ref_h_idx['#chr']]
                    ref_pos = int(ref_line[ref_h_idx['pos']])
                except IndexError:
                    ref_has_lines = False
            while ref_has_lines and int(ref_chr) == int(var.chr) and ref_pos == var.pos:
                ref_vars.append(Variant(ref_chr, ref_pos,
                                        ref_line[ref_h_idx['ref']], ref_line[ref_h_idx['alt']],
                                        ref_line[ref_h_idx['af_alt']]))
                ref_line = fp_ref.readline().strip().split('\t')
                try:
                    ref_chr = ref_line[ref_h_idx['#chr']]
                    ref_pos = int(ref_line[ref_h_idx['pos']])
                except IndexError:
                    ref_has_lines = False

            gnomad_af = 'NA'
            for r in ref_vars:
                if var.equalize_to(r):
                    gnomad_af = r.af
                    break

            print('\t'.join([var.chr, str(var.pos), var.ref, var.alt,
                             str(var.af), str(var.beta), str(var.se), str(var.pval),
                             str(gnomad_af)]))
            
def run():
    parser = argparse.ArgumentParser(description="Harmonize GWAS summary stats to reference")
    parser.add_argument('file_in', action='store', type=str, help='GWAS summary stats in SAIGE format')
    parser.add_argument('file_ref', action='store', type=str, help='Reference file')
    args = parser.parse_args()
    harmonize(args.file_in, args.file_ref)
              
if __name__ == '__main__':
    run()
