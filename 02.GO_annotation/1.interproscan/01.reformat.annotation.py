### This script is to combine multiple GO terms for each gene
## Chongjing Xia, 20240306


import pandas as pd
import sys

def reformat_file(input_file, output_file):
    # Load the data into a pandas DataFrame
    df = pd.read_csv(input_file, sep='\t', header=None, names=['gene', 'ipr', 'desc', 'go'])

    # Group by 'gene' and aggregate the other columns
    df_grouped = df.groupby('gene').agg({
        'ipr': lambda x: ';'.join(x),
        'desc': lambda x: ';'.join(x),
        'go': lambda x: '|'.join(x.dropna())
    })

    # Write the output to a file
    df_grouped.to_csv(output_file, sep='\t')

# Get the input and output file paths from command line arguments
input_file = sys.argv[1]
output_file = sys.argv[2]

# Call the function with your input and output file paths
reformat_file(input_file, output_file)

