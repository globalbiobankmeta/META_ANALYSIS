#!/usr/bin/env python3

import argparse
import json
import re

def run():

    parser = argparse.ArgumentParser(description='Create a json config file for meta-analysis')
    parser.add_argument('in_files_loc', action='store', type=str, help='A tab-delimited text file with gs:// locations of sumstat input files. One row is one meta')
    parser.add_argument('in_phenos_loc', action='store', type=str, help='A text file with phenotype names according to in_files_loc')
    args = parser.parse_args()

    with open(args.in_phenos_loc) as f:
        phenos = [line.strip().split('\t')[0] for line in f.readlines()]

    with open(args.in_files_loc) as f:
        if len(f.readlines()) != len(phenos):
            raise Exception("Different number of lines in " + args.in_phenos_loc + " and " + args.in_files_loc)

    with open(args.in_files_loc) as f:
        index = 0
        for line in f:
            files = line.strip().split('\t')
            conf = {'meta': []}
            conf_loo = [{'meta': []} for file in files]
            conf_leave_ukb = {'meta': []}
            n_cases_total = 0
            n_controls_total = 0
            print(phenos[index])
            for i, file in enumerate(files):
                filename = file.split('/').pop()
                cleanname = re.sub('\.munged(.)*$', '', filename)
                cleanname = re.sub('.formatted|.txt|.gz|.bgz', '', cleanname)
                fields = cleanname.split('.')
                if (len(fields) != 10):
                    raise Exception('Unexpected filename in {}: {}'.format(args.in_files_loc, filename))
                try:
                    n_cases = int(fields[6])
                    n_cases_total = n_cases_total + n_cases
                    n_controls = int(fields[7])
                    n_controls_total = n_controls_total + n_controls
                except ValueError:
                    raise Exception('Could not parse numbers of cases and controls from {}'.format(filename))
                cohort = {
                    'name': fields[0] + '_' + fields[5],
                    'file': file.replace('gs://', '/cromwell_root/'),
                    'n_cases': n_cases,
                    'n_controls': n_controls,
                    'chr':'#CHR',
                    'pos':'POS',
                    'ref':'Allele1',
                    'alt':'Allele2',
                    'effect':'BETA',
                    'effect_type':'beta',
                    'pval':'p.value',
                    'se':'SE',
                    'extra_cols':['AF_Allele2', 'AF_fc', 'imputationInfo', 'N']
                }
                conf['meta'].append(cohort)
                if not cohort['name'].startswith('UKBB'):
                    conf_leave_ukb['meta'].append(cohort)
                for ci, c in enumerate(conf_loo):
                    if ci != i:
                        c['meta'].append(cohort)
                    else:
                        c['cohort'] = cohort['name']
                print('{}_{}\t{}\t{}'.format(fields[0], fields[5], n_cases, n_controls))
            with open(phenos[index] + '.json', 'w') as out:
                json.dump(conf, out, indent=4)
            if len(conf['meta']) > len(conf_leave_ukb['meta']):
                with open(phenos[index] + '_leave_UKBB.json', 'w') as out:
                    json.dump(conf_leave_ukb, out, indent=4)
            for c in conf_loo:
                with open(phenos[index] + '_leave_' + c['cohort'] + '.json', 'w') as out:
                    json.dump({'meta': c['meta']}, out, indent=4)
            index = index + 1
            print()

if __name__ == '__main__':
    run()
