#! /usr/bin/env python
"""
Dig into shared nt segments by doing nucmer alignments.
"""
import argparse
import glob
import csv
import os
from pymummer import coords_file, alignment, nucmer
import gzip
from collections import defaultdict
import screed
import shutil
import math


# should change this so it's provided on command line, maybe...
GENOME_SUFFIXES = ['', '.fna', '.fa', '*_genomic.fna.gz', '*.fna', '*.fa']


def find_genome_filename(genomes_dir, ident):
    """Find the genome for the given identifier by guessing at extensions.
    Complain bitterly if a single match cannot be found."""

    for suffix in GENOME_SUFFIXES:
        pattern = os.path.join(genomes_dir, ident + suffix)
        matches = glob.glob(pattern)
        if len(matches) > 1:
            assert 0, "more than one match to {}; {}".format(ident, matches)

        if matches:
            assert len(matches) == 1
            return matches[0]

    print('no matches in {} for {}'.format(genomes_dir, ident))
    print('looking for suffixes:', GENOME_SUFFIXES)
    assert 0, "cannot find matches to {} in {}".format(ident, genomes_dir)


def remove_contigs(ident, genomefile, keep_d, verbose=True):
    """
    remove contigs from 'genomefile' whose names are in keep_d keys.

    removed contigs go in genomefile + '.removed.fa'
    retained contigs go in genomefile + '.kept.fa'

    return the total number of bp removed.
    """
    bp_skipped = 0
    contigs_skipped = 0
    bp_total = 0
    contigs_total = 0

    kept_outfp = open(genomefile + '.kept.fa', 'wt')
    removed_outfp = open(genomefile + '.removed.fa', 'wt')
    
    for record in screed.open(genomefile):
        bp_total += len(record.sequence)
        contigs_total += 1

        name = record.name.split(' ')[0]
        if name in keep_d:
            # filter
            bp_skipped += len(record.sequence)
            contigs_skipped += 1
            removed_outfp.write('>{}\n{}\n'.format(record.name, record.sequence))
        else:
            kept_outfp.write('>{}\n{}\n'.format(record.name, record.sequence))

    if verbose:
        print('{}: removed {:.0f}kb of {:.0f}kb ({:.0f}%), {} of {} contigs'.format(ident, bp_skipped / 1000, bp_total / 1000, bp_skipped / bp_total * 100, contigs_skipped, contigs_total))

    return bp_skipped


def copy_and_gunzip_genome(from_name, to_name):
        xopen = open
        if from_name.endswith('.gz'):
            xopen = gzip.open
        with xopen(from_name, 'rb') as fp1:
            with open(to_name, 'wb') as fp2:
                      fp2.write(fp1.read())


def main():
    p = argparse.ArgumentParser()
    p.add_argument('oddities_csv')
    p.add_argument('genomes_dir', help='fastani database dir')
    p.add_argument('--percent-threshold', type=float,
                   default=95.0)
    p.add_argument('--length-threshold', type=int, default=0)
    p.add_argument('-v', '--verbose', action='store_true')
    args = p.parse_args()

    print('loading', args.oddities_csv)
    print('getting genomes from:', args.genomes_dir)
    print('length threshold for alignments (bp):', args.length_threshold)
    print('lower cutoff for identity (%):', args.percent_threshold)

    prefix = args.oddities_csv
    assert prefix.endswith('.csv')
    if prefix.endswith('.csv'):
        prefix = prefix[:-4]

    fp = open(args.oddities_csv, 'rt')
    r = csv.DictReader(fp)

    alignments_dir = prefix + '.alignments'
    print('putting alignments in:', alignments_dir)
    try:
        os.mkdir(alignments_dir)
    except FileExistsError:
        print('warning: directory already exists!')

    print('----')

    for row in r:
        cluster_name = row['cluster']
        ident1 = os.path.basename(row['ident1'])
        ident2 = os.path.basename(row['ident2'])

        if args.verbose:
            print(cluster_name, ident1, ident2)

        # copy & name genome files "clusterx.y.IDENT.fa. gunzip if necessary,
        # since nucmer doesn't handle gzip.
        fn1 = find_genome_filename(args.genomes_dir, ident1)
        genome1 = os.path.join(alignments_dir, '{}.{}.fa'.format(cluster_name, ident1))
        copy_and_gunzip_genome(fn1, genome1)

        fn2 = find_genome_filename(args.genomes_dir, ident2)
        genome2 = os.path.join(alignments_dir, '{}.{}.fa'.format(cluster_name, ident2))
        copy_and_gunzip_genome(fn2, genome2)

        nucmer_output_name = os.path.join(alignments_dir, cluster_name + '.a')

        if not os.path.exists(nucmer_output_name):
            print('running {} alignments...'.format(cluster_name))
            runner = nucmer.Runner(genome1, genome2, nucmer_output_name)
            runner.run()
            print('...done!')
        else:
            if args.verbose:
                print('using cached alignments file', nucmer_output_name)

        file_reader = coords_file.reader(nucmer_output_name)
        alignments = [coord for coord in file_reader if not coord.is_self_hit()]

        # alignment obj:
        # 'frame', 'hit_length_qry', 'hit_length_ref', 'intersects_variant', 'is_self_hit', 'on_same_strand', 'percent_identity', 'qry_coords', 'qry_coords_from_ref_coord', 'qry_end', 'qry_length', 'qry_name', 'qry_start', 'ref_coords', 'ref_coords_from_qry_coord', 'ref_end', 'ref_length', 'ref_name', 'ref_start', 'reverse_query', 'reverse_reference', 'to_msp_crunch']

        # sort alignments by length of hit
        alignments.sort(key = lambda x: -x.hit_length_qry)

        # track alignments over a particular threshold
        keep_alignments = []
        all_bp = 0
        aligned_bp = 0
        weighted_percent_identity = 0.
        skipped_bp = 0
        skipped_aln = 0

        for alignment in alignments:
            weighted_percent_identity += alignment.percent_identity * alignment.hit_length_qry
            all_bp += alignment.hit_length_qry

            # do we pass the length and percent identitiy thresholds? if so,
            # keep!
            if alignment.hit_length_qry >= args.length_threshold and \
               alignment.percent_identity >= args.percent_threshold:
                aligned_bp += alignment.hit_length_qry
                keep_alignments.append(alignment)
            else:
                skipped_bp += alignment.hit_length_qry
                skipped_aln += 1

        # ditch if no alignments
        if not keep_alignments:
            print('** FLAG: no kept alignments for {}, punting.'.format(cluster_name))
            print('')
            continue

        # set up the printed out info
        lca_name = "(root)"     # if empty lca, => root of taxonomy.
        if row['lca']:
            lca_name = row['lca']

        shared_kmers = int(row['shared_kmers'])
        ksize = int(row['ksize'])

        # nice output! with some flags.
        print('{}: {:.0f}kb aln ({:.0f}k {}-mers) across {}; longest contig: {:.0f} kb'.format(cluster_name, aligned_bp / 1000, shared_kmers / 1000, ksize, lca_name, keep_alignments[0].hit_length_qry / 1000))
        print('weighted percent identity across alignments: {:.1f}%'.format(weighted_percent_identity / all_bp))
        print('skipped {:.0f} kb of alignments in {} alignments (< {} bp or < {:.0f}% identity)'.format(skipped_bp / 1000, skipped_aln, args.length_threshold, args.percent_threshold))
        if abs(math.log(shared_kmers / aligned_bp) > 1):
            print('** FLAG, oddly too little or too many aligned bp vs k-mers')

        ### track & remove contigs from query genome (genome2)

        keep_d = defaultdict(set)
        for aln in keep_alignments:
            keep_d[aln.qry_name].add(aln)

        bp_removed = remove_contigs(ident2, genome2, keep_d)

        flag_2 = 0
        if bp_removed > 2.5*aligned_bp:
            flag_2 = 1

            # reset to rm kept, and removed is empty.
            os.unlink(genome2 + '.kept.fa')
            with open(genome2 + '.removed.fa', 'wt') as fp:
                pass

        ### track & remove contigs from ref genome (genome1)

        keep_d = defaultdict(set)
        for aln in keep_alignments:
            keep_d[aln.ref_name].add(aln)

        bp_removed = remove_contigs(ident1, genome1, keep_d)

        flag_1 = 0
        if bp_removed > 2.5*aligned_bp:
            flag_1 = 1
            # reset to rm kept, and removed is empty.
            os.unlink(genome1 + '.kept.fa')
            with open(genome1 + '.removed.fa', 'wt') as fp:
                pass

        # output summary of flags

        if flag_1 and flag_2:
            print('** FLAGFLAG, too much removed from both!')
        elif flag_1 and not flag_2:
            print('** FLAG, {} is probably contaminated (too much rm from {})'.format(ident2, ident1))
        elif flag_2 and not flag_1:
            print('** FLAG, {} is probably contaminated (too much rm from {})'.format(ident1, ident2))

        print('')

        # done with main loop!


if __name__ == '__main__':
    main()
