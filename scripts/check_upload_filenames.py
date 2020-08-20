#!/usr/bin/python

import sys, re

POPS = ['EUR', 'AFR', 'FIN', 'ARAB', 'SAS', 'EAS', 'OTH', 'CSA', 'AMR', 'NON-EUR', 'ALL']
ANAS = ['A1', 'A2', 'B1', 'B2', 'B3', 'C1', 'C2', 'D1']
TYPES = ['ALL', 'MALE', 'FEMALE']

with open(sys.argv[1], 'rt') as f:
    for line in f:
        ok = True
        line = line.strip()
        filename = line.split('/').pop()
        filename = re.sub('.formatted|.txt|.gz|.bgz', '', filename)
        #filename = re.sub('\.formatted(.)*$', '', filename)
        fields = filename.split('.')
        if (len(fields) != 10):
            print(str(len(fields)) + ' fields\t' + line)
            continue
        try:
            freeze = int(fields[3])
            n_cases = int(fields[6])
            n_controls = int(fields[7])
            date = int(fields[9])
        except ValueError:
            print('Number error\t{}'.format(line))
            ok = False
        if fields[5] not in POPS:
            print('Unexpected pop\t{}'.format(line))
            ok = False
        ana_ok = False
        for ana in ANAS:
            if ana in fields[2]:
                ana_ok = True
                break
        if not ana_ok:
            print('Unexpected ANA\t{}'.format(line))
            ok = False
        if fields[4].upper() not in TYPES:
            print('Unexpected type\t{}'.format(line))
            ok = False
        if ok:
            print('ALLGOOD\t{}'.format(line))
