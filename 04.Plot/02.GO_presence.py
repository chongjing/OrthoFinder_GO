#!/usr/bin/env python3
"""
Convert multiple GO term lists into a presence/absence matrix.

Usage:
    python3 format.py 01.AllGO.list DEGs_*_log2FC_up_GO.csv > output.csv

Arguments:
    01.AllGO.list              - File containing all possible GO terms (one per line)
    DEGs_*_log2FC_up_GO.csv    - One or more species-specific GO term files.
                                 Species name is extracted from the part between
                                 "DEGs_" and "_log2FC_up_GO.csv".

Output:
    CSV with first column = species name, subsequent columns = GO terms from the
    full list, each cell = 1 if the term is present in that species, 0 otherwise.
"""

import sys
import os

def extract_species_name(filename):
    """Extract species name from a filename like 'DEGs_mouse_log2FC_up_GO.csv'."""
    basename = os.path.basename(filename)
    # Remove prefix and suffix
    if basename.startswith("DEGs_") and basename.endswith("_log2FC_up_GO.csv"):
        species = basename[5:-19]  # len("DEGs_")=5, len("_log2FC_up_GO.csv")=19
        return species
    else:
        # If filename does not match expected pattern, use the whole name (without extension) as fallback
        return os.path.splitext(basename)[0]

def main():
    if len(sys.argv) < 3:
        print("Error: Insufficient arguments.", file=sys.stderr)
        print(__doc__, file=sys.stderr)
        sys.exit(1)

    full_list_file = sys.argv[1]
    species_files = sys.argv[2:]

    # Read full list of GO terms (preserve order)
    try:
        with open(full_list_file, 'r') as f:
            go_terms = [line.strip() for line in f if line.strip()]
    except IOError as e:
        print(f"Error reading full list file: {e}", file=sys.stderr)
        sys.exit(1)

    if not go_terms:
        print("Warning: Full GO list is empty.", file=sys.stderr)

    # Collect data for each species
    species_data = []  # list of tuples (species_name, set_of_terms)
    for sf in species_files:
        species = extract_species_name(sf)
        try:
            with open(sf, 'r') as f:
                # Read species-specific GO terms into a set for fast lookup
                terms = {line.strip() for line in f if line.strip()}
            species_data.append((species, terms))
        except IOError as e:
            print(f"Warning: Cannot read {sf}, skipping. Error: {e}", file=sys.stderr)
            continue

    # Output CSV
    # Header
    header = ["Species"] + go_terms
    print(",".join(header))

    # Data rows
    for species, terms_set in species_data:
        row = [species] + ["1" if go in terms_set else "0" for go in go_terms]
        print(",".join(row))

if __name__ == "__main__":
    main()
