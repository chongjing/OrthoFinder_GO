### this script is to extract orthologroup from re-formated orthofinder results and append it to DE results
#Chongjing Xia 20240306

import pandas as pd
import sys

def main(orthologous_file, gene_expression_file, output_file):
    # Read the orthologous group file into a pandas DataFrame
    df_ortho = pd.read_csv(orthologous_file, sep='\t')

    # Create a dictionary mapping gene names to orthologous groups
    gene_to_group = {}
    for _, row in df_ortho.iterrows():
        group = row[0]
        for genes in row[1:]:
            if isinstance(genes, str):
                for gene in genes.split(','):
                    gene_to_group[gene] = group

    # Read the gene expression file into a pandas DataFrame
    df_gene = pd.read_csv(gene_expression_file, sep='\t')

    # Add a new column for the orthologous group
    df_gene['Orthologous Group'] = df_gene.iloc[:, 0].map(gene_to_group)

    # Write the DataFrame to a new file
    df_gene.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    # The first command line argument is the orthologous group file,
    # the second is the gene expression file, and the third is the output file
    main(sys.argv[1], sys.argv[2], sys.argv[3])
