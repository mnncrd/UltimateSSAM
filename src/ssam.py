"""Ultimate Secondary Structure Assignment Method.

This program assigns secondary structures to a sequence.

    Typical usage example:

    $ python ssam.py filename -o output_file
"""

import argparse

class Atom():

    """
    Summary of class Atom.

    Lorem ipsum.

    Attributes:
        code: A dictionnary mapping the 3-letter code to the 1-letter code
            for the Residue name.
        atom_name: A dictionary of Atoms instance forming the Residue.
        aa_name: A string indicating the name of the Residue.
        aa_chain: A string indicating on which chain is the Residue.
        aa_nb: A integer indicating the number of the Residue.
        x_coord: A float representing the x coordinate of the atom.
        y_coord: A float representing the y coordinate of the atom.
        z_coord: A float representing the z coordinate of the atom.
    """

    code = {
        "ALA":"A", "CYS":"C", "ASP":"D", "GLU":"E", "PHE":"F", "GLY":"G",
        "HIS":"H", "ILE":"I", "LYS":"K", "LEU":"L", "MET":"M", "ASN":"N",
        "PRO":"P", "GLN":"Q", "ARG":"R", "SER":"S", "THR":"T", "VAL":"V",
        "TRP":"W", "TYR":"Y"
    }

    def __init__(self, line):

        self.atom_name = line[12:16].strip()
        self.aa_name = self.code[line[17:20]]
        self.aa_chain = line[21:22]
        self.aa_nb = int(line[22:26])
        self.x_coord = float(line[30:38])
        self.y_coord = float(line[38:46])
        self.z_coord = float(line[46:54])

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
