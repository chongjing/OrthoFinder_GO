### This script is to label genes within each orthologroup "UP","DOWN", or "NON" expressed
### Chongjing Xia (cx264@cam.ac.uk), 20240502
### Usage: python 00.ortho_express_summery.py DEGs_arab_log2FC_up DEGs_arab_log2FC_down 02.Orthogroups.Arab.Express_label.tsv 

import csv
import sys

# Function to read genes and their expression status from a file
def read_genes(filename, status):
    genes = {}
    with open(filename, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            gene = row['Gene'] if 'Gene' in row else row['Genes']
            genes[gene] = status
    return genes

# Main function to process the files and generate the output
def process_files(orthogroups_file, up_file, down_file, output_file):
    # Read up- and down-expressed genes
    up_genes = read_genes(up_file, 'UP')
    down_genes = read_genes(down_file, 'DOWN')

    # Open the orthogroups file and the output file
    with open(orthogroups_file, 'r') as ortho_file, open(output_file, 'w', newline='') as out_file:
        ortho_reader = csv.reader(ortho_file, delimiter='\t')
        out_writer = csv.writer(out_file, delimiter='\t')

        # Write the header row to the output file
        header = next(ortho_reader)
        out_writer.writerow(header)

        # Process each orthogroup starting from the second row
        for row in ortho_reader:
            orthogroup = row[0]
            if len(row) > 1 and row[1]:
                # Split the genes and determine their expression status
                genes = row[1].split(',')
                labeled_genes = []
                for gene in genes:
                    if gene in up_genes:
                        labeled_genes.append(f"{gene}(UP)")
                    elif gene in down_genes:
                        labeled_genes.append(f"{gene}(DOWN)")
                    else:
                        labeled_genes.append(f"{gene}(NON)")
                # Write the updated orthogroup with labeled genes to the output
                out_writer.writerow([orthogroup, ','.join(labeled_genes)])
            else:
                # Write orthogroups with no genes unchanged
                out_writer.writerow(row)

# Check if the script is run with the correct number of arguments
if len(sys.argv) != 5:
    print("Usage: python express_summerize.py <01.Orthogroups.Arab.tsv> <DEGs_arab_log2FC_up> <DEGs_arab_log2FC_down> <output.tsv>")
else:
    process_files(*sys.argv[1:])

