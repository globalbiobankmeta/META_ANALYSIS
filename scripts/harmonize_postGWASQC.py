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

def format_num(num, precision=4):
    return numpy.format_float_scientific(num, precision=precision) if num is not None else "NA"

class Variant():

    def __init__(self, chr, pos, ref, alt, af, filt, an):
        self.chr = chr
        self.pos = int(float(pos))
        self.ref = ref.strip().upper()
        self.alt = alt.strip().upper()
        self.af = float(af) if af != "NA" else None
        self.filt = filt
        self.an = int(an)

    def __eq__(self, other):

        return self.chr == other.chr and self.pos == other.pos and self.ref == other.ref and self.alt == other.alt

    def __lt__(self, other):

        return (  (self.chr==other.chr and self.pos<other.pos)
                  or (self.chr < other.chr)
               )

    def is_equal(self, other:'Variant') -> bool:
        """
            Checks if this Variant is the same variant (possibly different strand or ordering of alleles)
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

    def __init__(self, chr, pos, ref, alt, extra_cols=[]):
        self.chr = chr
        self.pos = int(float(pos))
        self.ref = ref.strip().upper()
        self.alt = alt.strip().upper()
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
                    #self.beta = -1 * self.beta if self.beta is not None else None
                    #self.af = 1 - self.af if self.af is not None else None
                    #t = self.alt
                    #self.alt = self.ref
                    #self.ref = t
                    return True

            elif (self.ref == other.alt and self.alt == other.ref) :
                #self.beta = -1 * self.beta if self.beta is not None else None
                #self.af = 1 - self.af if self.af is not None else None
                #t = self.alt
                #self.alt = self.ref
                #self.ref = t
                return True
            elif (self.ref == flip_ref and self.alt==flip_alt):
                #self.ref = flip_strand(self.ref)
                #self.alt = flip_strand(self.alt)
                return True
            elif (self.ref == flip_alt and self.alt==flip_ref):
                #self.beta = -1 * self.beta if self.beta is not None else None
                #self.af = 1 - self.af if self.af is not None else None
                #self.ref =flip_strand(self.alt)
                #self.alt = flip_strand(self.ref)
                return True

        return False

def harmonize(file_in, file_remove, file_ambiguity):

    
    rmd = {}
    with open(file_remove, "r") as fp_rm:
        headtemp = fp_rm.readline()
        for line in fp_rm.readlines():
            if 'this' in line: 
                linea = line.strip().split()
                key = linea[0]+":"+linea[1]+":"+linea[2]+":"+linea[3]
                rmd[key] = key


    
    amd = {}
    with open(file_ambiguity, "r") as fp_am:
        headtemp = fp_am.readline()
        for line in fp_am.readlines():
            if 'this' in line: 
                linea = line.strip().split()
                key = linea[0]+":"+linea[1]+":"+linea[2]+":"+linea[3]
                amd[key] = key



    with gzip.open(file_in, 'rt') as f:

        headerline0 = f.readline().strip()
        headerline = headerline0.split('\t')
        h_idx = {h:i for i,h in enumerate(headerline)}
	
        #print('#CHR\tPOS\tAllele1\tAllele2\tAF_Allele2\timputationInfo\tBETA\tSE\tp.value\tAF_' + file_ref.replace('.gz', '') + '\tAF_fc\tN')
        print(headerline0 + "\tstrandflip")
        for line in f:
            s = line.strip().split('\t')
            var = VariantData(s[h_idx['#CHR']].replace('chr', '').replace('X', '23'), s[h_idx['POS']],
                              s[h_idx['Allele1']], s[h_idx['Allele2']])
            keyval = str(var.chr)+":"+str(var.pos)+":"+ var.ref+":"+var.alt
            if keyval not in rmd: 
                if keyval in amd:   
                    s[4] = str(1-float(s[4]))
                    s[6] = str(-1*float(s[6]))
                    print('\t'.join(s) +"\tTRUE")
                else:
                    print('\t'.join(s) +"\tFALSE")
                     
            
def run():
    parser = argparse.ArgumentParser(description="Harmonize GWAS summary stats to reference")
    parser.add_argument('file_in', action='store', type=str, help='GWAS summary stats in SAIGE format')
    parser.add_argument('file_remove', action='store', type=str, help='remove file')
    parser.add_argument('file_ambiguity', action='store', type=str, help='ambuguity file')
    args = parser.parse_args()
    harmonize(args.file_in, args.file_remove, args.file_ambiguity)
    
if __name__ == '__main__':
    run()
