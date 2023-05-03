## python3
import numpy as np
import pandas as pd
import sys
import argparse

parser = argparse.ArgumentParser(description='Filter a phylip file by individual-level missing data',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-I', '--infile', type=str, help='input file path')
parser.add_argument('-O', '--outfile', type=str, help='output file path')
parser.add_argument('-M', '--miss', type=float, help='specify the minimum proportion of missing variants per individual')

args = vars(parser.parse_args())

infile=args['infile']
outfile=args['outfile']
ind_missing_threshold=args['miss']

# Read in the phylip file

data = pd.read_csv(infile, delim_whitespace=True, skiprows=1, header=None) # read in dataframe, skip header
data.columns = ['ID','seq']

phylip_header = pd.read_csv(infile,  delim_whitespace=True, header=None, skiprows = lambda x: x not in [0]) # take the phylip header; heavy-handed but it wokrs
phylip_header.columns = ['sequences', 'length']

# Ensure the the header is correct
if len(data)!=phylip_header['sequences'][0]:
  print("infile formatted incorrectly")

if len(data['seq'][0])!=phylip_header['length'][0]:
  print("infile formatted incorrectly")


length=phylip_header['length'][0] # extract sequence length form phylip header
IDs_to_remove=list() # empty list to fill

# Look through the sequences, marking IDs wtih missing data above the threshold
for i,row in data.iterrows():
    ind_missing_count_i=data['seq'][i].count('N')
    ind_missing_proportion_i=ind_missing_count_i/length
    if ind_missing_proportion_i>ind_missing_threshold:
        IDs_to_remove.append(data['ID'][i])
    else:
        pass

# Remove IDs that breach the threshold
data_filtered = pd.DataFrame(data[~data['ID'].isin(IDs_to_remove)])

print(len(data_filtered),'of',phylip_header['sequences'][0],'sequences remain after filtering')

# Update the phylip header to reflect the new number of individuals 
new_phylip_header = pd.DataFrame(columns = ['ID', 'seq'], index=[0])
new_phylip_header['ID'][0]=len(data_filtered)
new_phylip_header['seq'][0]=phylip_header['length'][0]
output = pd.concat([new_phylip_header, data_filtered], ignore_index=True)

output.to_csv(outfile, sep=' ', index=False, header=False) # write output
