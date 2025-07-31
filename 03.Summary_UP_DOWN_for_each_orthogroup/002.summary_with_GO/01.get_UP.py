import sys

def filter_genes(input_file, output_file):
    """
    Reads the input file, filters out genes with labels 'DOWN' and 'NON',
    and writes the result to the output file.

    Args:
    input_file (str): Path to the input file.
    output_file (str): Path to the output file.
    """
    with open(input_file, 'r') as infile:
        with open(output_file, 'w') as outfile:
            for line in infile:
                # Split the line into orthogroup name and genes list
                parts = line.strip().split('\t')
                orthogroup = parts[0]
                if len(parts) > 1:
                    genes = parts[1].split(',')

                    # Filter out genes with labels 'DOWN' and 'NON'
                    filtered_genes = [gene for gene in genes if '(UP)' in gene]

                    # Write the orthogroup and filtered genes to the output file
                    if filtered_genes:
                        outfile.write(f"{orthogroup}\t{','.join(filtered_genes)}\n")
                    else:
                        outfile.write(f"{orthogroup}\n")
                else:
                    # If no genes are listed, just write the orthogroup name
                    outfile.write(f"{orthogroup}\n")

if __name__ == "__main__":
    # Check if the correct number of command-line arguments were provided
    if len(sys.argv) != 3:
        print("Usage: python mycode.py input output")
        sys.exit(1)

    # Get the input and output file paths from the command-line arguments
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # Run the filter function
    filter_genes(input_file, output_file)
