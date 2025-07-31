import pandas as pd
import re
import sys

def process_genes(genes):
    """
    This function takes a string of comma-separated gene names,
    removes version numbers and duplicates from each gene name,
    and returns a string of the processed gene names.
    """
    if not isinstance(genes, str):
        return genes

    # Split the string into a list of gene names
    genes_list = genes.split(', ')
    
    # Use a regular expression to remove version numbers from each gene name
    processed_genes = [re.sub(r'\.\d+', '', gene) for gene in genes_list]
    
    # Remove duplicates by converting the list to a set, then convert it back to a list
    processed_genes = list(set(processed_genes))
    
    # Join the processed gene names back into a string
    return ','.join(processed_genes)

def main(input_file, output_file):
    # Read the TSV file into a pandas DataFrame
    df = pd.read_csv(input_file, sep='\t')

    # Apply the process_genes function to each cell in the DataFrame
    df = df.applymap(process_genes)

    # Write the processed DataFrame to a new TSV file
    df.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    # The first command line argument is the input file, and the second is the output file
    main(sys.argv[1], sys.argv[2])

