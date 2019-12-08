"""Ultimate Secondary Structure Assignment Method.

This program assigns secondary structures to a sequence.

    Typical usage example:

    $ python ssam.py filename -o output_file
"""

import argparse

def open_pdb_file(filename):

    """Reads the .pdb file and stores info"""

    assert filename.lower().endswith(".pdb"), "Program can only work with a .pdb file"
    print("Reading file")
    file_pdb = open(filename, "r")
    file_pdb.close()
    pdb_info = []
    residues = []
    print("ok")
    return pdb_info, residues

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
