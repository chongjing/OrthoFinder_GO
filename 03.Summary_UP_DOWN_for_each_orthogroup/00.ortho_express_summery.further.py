import csv
import sys

# Main function to summarize the expression labels for each orthogroup
def summarize_expression(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')

        # Write the header to the output file
        header = next(reader)
        writer.writerow(header)

        # Process each orthogroup to count the UP, DOWN, and NON labels
        for row in reader:
            orthogroup = row[0]
            # Initialize counters
            counts = {'UP': 0, 'DOWN': 0, 'NON': 0}

            if len(row) > 1 and row[1]:
                # Split the genes and count the expression status
                genes = row[1].split(',')
                for gene in genes:
                    label = gene.split('(')[-1].rstrip(')')
                    counts[label] += 1

                # Format the summary string
                summary = f"UP({counts['UP']});DOWN({counts['DOWN']});NON({counts['NON']})"
                writer.writerow([orthogroup, summary])
            else:
                # For orthogroups with no genes, write the row unchanged
                writer.writerow(row)

# Check if the script is run with the correct number of arguments
if len(sys.argv) != 3:
    print("Usage: python summarize_expression.py <input.tsv> <output.tsv>")
else:
    summarize_expression(*sys.argv[1:])

