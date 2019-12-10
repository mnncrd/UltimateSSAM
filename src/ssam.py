"""Ultimate Secondary Structure Assignment Method.

This program assigns secondary structures to a sequence.

    Typical usage example:

    $ python ssam.py filename -o output_file
"""
#import
import argparse
import math
#import modules
import angles
import vectors

class Residue():
    """
    Summary of class Residue.

    Lorem ipsum.

    Attributes:
        atoms: A dictionary of Atoms instance forming the Residue.
        number: A integer indicating the number of the Residue.
        name:  A string indicating the name of the Residue.
        chain: A string indicating on which chain is the Residue.
        angles: A dictionnary of angles values for the Residue.
    """

    def __init__(self, atoms, res_nb=None):

        self.atoms = {
            "CA":next((a for a in atoms if a.atom_name == "CA"), None),
            "C":next((a for a in atoms if a.atom_name == "C"), None),
            "H":next((a for a in atoms if a.atom_name == "H"), None),
            "O":next((a for a in atoms if a.atom_name == "O"), None),
            "N":next((a for a in atoms if a.atom_name == "N"), None)
        }
        if res_nb is None:
            self.number = self.atoms["CA"].aa_nb
        else:
            self.number = res_nb
        self.name = self.atoms["CA"].aa_name
        self.chain = self.atoms["CA"].aa_chain
        self.angles = {
            "TCO":"0.000",
            "KAPPA":"360.0",
            "ALPHA":"360.0",
            "PHI":"360.0",
            "PSI":"360.0"
        }
        self.struct = {
            "H":False, "B":False, "E":False, "G":False, "I":False, "T":False, "S":False,
            "3":" ", "4":" ", "5":" ",
            "BEND":" ",
            "CHR":" "
        }

    def compute_tco(self, oth):

        """Computes the tco"""

        v_i = self.atoms["C"].vector(self.atoms["O"])
        v_j = oth.atoms["C"].vector(oth.atoms["O"])
        angle = math.cos(angles.compute_angle(v_i, v_j))
        self.angles["TCO"] = "{:.3f}".format(angle)

    def compute_kappa(self, oth_b, oth_a):

        """Computes kappa"""

        v_i = self.atoms["CA"].pos_vector()
        v_j = oth_b.atoms["CA"].pos_vector()
        v_k = oth_a.atoms["CA"].pos_vector()
        v_ij = vectors.compute_diff_vect(v_i, v_j)
        v_ki = vectors.compute_diff_vect(v_k, v_i)
        angle = math.degrees(angles.compute_angle(v_ki, v_ij))
        self.angles["KAPPA"] = "{:.1f}".format(angle)
        if angle > 70:
            self.struct["S"] = True
            self.struct["BEND"] = "S"

    def compute_alpha(self, oth_b, oth_a, oth_aa):

        """Computes alpha"""

        v_i = oth_b.atoms["CA"].pos_vector()
        v_j = self.atoms["CA"].pos_vector()
        v_k = oth_a.atoms["CA"].pos_vector()
        v_l = oth_aa.atoms["CA"].pos_vector()
        angle = angles.compute_dihedral_angle(v_i, v_j, v_k, v_l)
        self.angles["ALPHA"] = "{:.1f}".format(angle)
        if 0 < angle < 180:
            self.struct["CHR"] = "+"
        elif -180 < angle < 0:
            self.struct["CHR"] = "-"

    def compute_phi(self, oth):

        """Computes phi"""

        v_i = oth.atoms["C"].pos_vector()
        v_j = self.atoms["N"].pos_vector()
        v_k = self.atoms["CA"].pos_vector()
        v_l = self.atoms["C"].pos_vector()
        angle = angles.compute_dihedral_angle(v_i, v_j, v_k, v_l)
        self.angles["PHI"] = "{:.1f}".format(angle)

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

    def compute_distance(self, oth):

        """Computes the distance between two atoms"""

        x_square = (self.x_coord - oth.x_coord)**2
        y_square = (self.y_coord - oth.y_coord)**2
        z_square = (self.z_coord - oth.z_coord)**2

        return math.sqrt(x_square + y_square + z_square)

    def vector(self, oth):

        """Returns the vector between two atoms"""

        vect_x = self.x_coord - oth.x_coord
        vect_y = self.y_coord - oth.y_coord
        vect_z = self.z_coord - oth.z_coord
        coords = [vect_x, vect_y, vect_z]
        return coords

    def pos_vector(self):

        """Returns the position vector of the atom"""

        vect_x = self.x_coord
        vect_y = self.y_coord
        vect_z = self.z_coord
        coords = [vect_x, vect_y, vect_z]
        return coords

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
        elif line[0:4] == 'ATOM':
            atom = Atom(line)
            if atom.aa_nb != res_nb:
                if len(atoms) > 0:
                    res = Residue(atoms)
                    residues.append(res)
                    atoms = []
                res_nb = atom.aa_nb
            atoms.append(atom)
    res = Residue(atoms)
    residues.append(res)
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

if __name__ == '__main__':
    main()
