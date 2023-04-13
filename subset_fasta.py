#! /usr/bin/env python3
import os
from Bio import SeqIO
import argparse
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import warnings
warnings.simplefilter("ignore")


parser = argparse.ArgumentParser(description='Subset a region of interest from a fasta file',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-I', '--indir', type=str, help='path to input directory')
parser.add_argument('-O', '--outdir', type=str, help='path to output directory')
parser.add_argument('-F', '--fasta', type=str, help='input fasta file')
parser.add_argument('-S', '--start', type=int, help='start position of region of interest (1-indexed)')
parser.add_argument('-E', '--end', type=int, help='end position of region of interest (1-indexed)')
parser.add_argument('-N', '--name', type=str, help='name (will become ID in output fasta file)')

args = vars(parser.parse_args())

# Subset a region of interest (e.g. a gene) from a larger .fa file
INDIR=args['indir']
OUTDIR=args['outdir']
fasta=args['fasta']
start=args['start']
end=args['end']
name=args['name']

for scaffold in SeqIO.parse(os.path.join(INDIR, fasta), 'fasta'):
    seq=scaffold.seq

region_of_interest=str(seq[start-1:end-1]) # Python is zero-indexed, so subtract 1 from the start and end positions
region_ID = ('>' + name) # 
fa_components = [region_ID, region_of_interest]

outfile=os.path.join(OUTDIR, ('{}.fa'.format(name)))

with open(outfile, 'w') as file:
    file.write('\n'.join(fa_components))