import sys

def parse_file1(file1):
    """Parses the first file and returns a dictionary with orthogroup ID as keys and lists of gene IDs as values."""
    orthogroups = {}
    with open(file1, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            orthogroup_id = parts[0]
            if len(parts) > 1 and parts[1]:
                genes = [gene.split('(')[0] for gene in parts[1].split(',')]
            else:
                genes = []
            orthogroups[orthogroup_id] = genes
    return orthogroups

def parse_file2(file2):
    """Parses the second file and returns a dictionary with gene IDs as keys and GO terms as values."""
    gene_go_terms = {}
    with open(file2, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
#            gene_id = parts[0].split('.')[0]  # Remove transcript version
            gene_id = parts[0] # No need to split on '.' since gene IDs match directly
            go_terms = parts[1] if len(parts) > 1 else ''
            gene_go_terms[gene_id] = go_terms
    return gene_go_terms

def map_go_terms(orthogroups, gene_go_terms):
    """Maps GO terms from gene_go_terms to orthogroups and returns a dictionary with orthogroup ID as keys and GO terms as values."""
    orthogroup_go_terms = {}
    for orthogroup_id, genes in orthogroups.items():
        go_terms_set = set()
        for gene in genes:
            if gene in gene_go_terms:
                go_terms_set.update(gene_go_terms[gene].split(','))
        orthogroup_go_terms[orthogroup_id] = ','.join(sorted(go_terms_set)) if go_terms_set else ''
    return orthogroup_go_terms

def write_output_file(out_file, orthogroup_go_terms):
    """Writes the output to a file with orthogroup ID and corresponding GO terms."""
    with open(out_file, 'w') as f:
        for orthogroup_id, go_terms in orthogroup_go_terms.items():
            f.write(f"{orthogroup_id}\t{go_terms}\n")

def main(file1, file2, out_file):
    orthogroups = parse_file1(file1)
    gene_go_terms = parse_file2(file2)
    orthogroup_go_terms = map_go_terms(orthogroups, gene_go_terms)
    write_output_file(out_file, orthogroup_go_terms)

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python mycode.py first.file second.file out.file")
        sys.exit(1)

    file1 = sys.argv[1]
    file2 = sys.argv[2]
    out_file = sys.argv[3]

    main(file1, file2, out_file)
