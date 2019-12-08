"""Ultimate Secondary Structure Assignment Method.

This program assigns secondary structures to a sequence.

    Typical usage example:

    $ python ssam.py filename -o output_file
"""

import argparse

def read_pdb_file(lines):

    """Stores info"""

    header_pdb = lines[0].strip()
    authors = []
    organism = ""
    molecule = ""
    residues = []
    res_nb = 0
    atoms = []
    for line in lines:
        if 'ORGANISM_SCIENTIFIC: ' in line:
            organism += line.strip()
        elif 'MOLECULE: ' in line:
            molecule += line.strip()
        elif line[0:6] == 'AUTHOR':
            authors.append(line[10:].strip().split(","))
    authors = ','.join(sum(authors, []))
    pdb_info = (header_pdb, organism, molecule, authors)
    return pdb_info, residues

def open_pdb_file(filename):

    """Opens the .pdb file and read the lines"""

    assert filename.lower().endswith(".pdb"), "Program can only work with a .pdb file"
    try:
        with open(filename, "r") as file_pdb:
            lines = file_pdb.readlines()
            pdb_info, residues = read_pdb_file(lines)
        file_pdb.close()
        print("ok")
        return pdb_info, residues
    except FileNotFoundError as fnf_error:
        print(fnf_error)

def main():

    """Main function"""

    #Get the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", help="The file to assign a secondary structure to", type=str)
    parser.add_argument("-o", type=str)
    parser.add_argument("-hy", "--hydrogen", action="store_true", help="add hydrogen atoms")
    args = parser.parse_args()

    #Read the file
    try:
        pdb_info, residues = open_pdb_file(args.filename)
    except AssertionError as error:
        print(error)
        print(args.filename, " is not a .pdb file")
