#! /usr/bin/env python
import os
import csv
import sys
import argparse
import pprint

import sourmash
from sourmash.lca import lca_utils


def main():
    p = argparse.ArgumentParser()
    p.add_argument('gtdbtk_dir')
    p.add_argument('lineages_csv_out')
    p.add_argument('-v', '--verbose', action='store_true')
    p.add_argument('--filter-prefix', default='')
    args = p.parse_args()

    tk_d = {}
    bac_summary = os.path.join(args.gtdbtk_dir, 'gtdbtk.bac120.summary.tsv')
    if os.path.exists(bac_summary):
        with open(bac_summary, 'rt') as fp:
            r = csv.DictReader(fp, delimiter='\t')
            for row in r:
                tk_d[row['user_genome']] = row['classification']

    ar_summary = os.path.join(args.gtdbtk_dir, 'gtdbtk.ar122.summary.tsv')
    if os.path.exists(ar_summary):
        with open(ar_summary, 'rt') as fp:
            r = csv.DictReader(fp, delimiter='\t')
            for row in r:
                tk_d[row['user_genome']] = row['classification']

    print('loaded {} rows from gtdbtk classify_wf in dir {}'.format(len(tk_d), args.gtdbtk_dir))

    lca_d = {}
    fp = open(args.lineages_csv_out, 'wt')
    w = csv.writer(fp)

    w.writerow('accession,superkingdom,phylum,class,order,family,genus,species'.split(','))

    n_skipped = 0
    n_written = 0
    for ident, tax in tk_d.items():
        if args.filter_prefix and not ident.startswith(args.filter_prefix):
            n_skipped += 1
            continue

        row = [ident] + tax.split(';')
        while len(row[-1]) == 3 and row[-1].endswith('__'):
            row.pop()
        w.writerow(row)
        n_written += 1

    if not n_written:
        raise Exception("0 written of {} rows - something went wrong!".format(len(tk_d)))

    print('wrote {} rows to {}'.format(n_written, args.lineages_csv_out))
    if n_skipped:
        print('(skipped {} because of --filter-prefix)'.format(n_skipped))


if __name__ == '__main__':
    main()
