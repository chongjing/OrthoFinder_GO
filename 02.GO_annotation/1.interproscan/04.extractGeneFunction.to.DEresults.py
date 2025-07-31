## this script is to extract gene function from formated Interproscan results and append it to DE results
#Chongjing Xia 20240306

import pandas as pd
import sys

def main(gene_file, function_file, output_file):
    # Read the gene file into a pandas DataFrame
    df_gene = pd.read_csv(gene_file, sep='\t')

    # Read the function file into a pandas DataFrame
    df_function = pd.read_csv(function_file, sep='\t')

    # Ensure the gene name column has the same name in both dataframes
    df_gene.rename(columns={df_gene.columns[0]: 'Gene'}, inplace=True)
    df_function.rename(columns={df_function.columns[0]: 'Gene'}, inplace=True)

    # Merge the two dataframes on the gene name column
    df_merged = pd.merge(df_gene, df_function, on='Gene', how='left')

    # Write the merged DataFrame to a new file
    df_merged.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    # The first command line argument is the gene file,
    # the second is the function file, and the third is the output file
    main(sys.argv[1], sys.argv[2], sys.argv[3])
