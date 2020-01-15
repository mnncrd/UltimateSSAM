"""Ultimate Secondary Structure Assignment Method.

This program assigns secondary structures to a sequence.

    Typical usage example:

    $ python ssam.py filename -o output_file
"""
# Import libraries
import argparse
import math
import sys

# Import modules
import angles
import vectors
import dssp
import checks

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

    def __init__(self, atoms):

        self.atoms = {
            "CA":next((a for a in atoms if a.atom_name == "CA"), None),
            "C":next((a for a in atoms if a.atom_name == "C"), None),
            "H":next((a for a in atoms if a.atom_name == "H"), None),
            "O":next((a for a in atoms if a.atom_name == "O"), None),
            "N":next((a for a in atoms if a.atom_name == "N"), None)
        }
        self.number = atoms[0].aa_nb
        self.name = atoms[0].aa_name
        self.chain = atoms[0].aa_chain
        self.angles = {
            "TCO":0.000,
            "KAPPA":360.0,
            "ALPHA":360.0,
            "PHI":360.0,
            "PSI":360.0
        }
        self.struct = {
            "STRC":" ",
            "H":False, "B":False, "E":False, "G":False, "I":False, "T":False, "S":False,
            "3":" ", "4":" ", "5":" ",
            "BEND":"",
            "CHR":"",
            "SHEET": "",
            "LAD1":"", "LAD2":"",
            "BP1":0, "BP2":0
        }
        self.bonds = {
            "O1": 0, "VO1": 0,
            "O2": 0, "VO2": 0,
            "N1": 0, "VN1": 0,
            "N2": 0, "VN2": 0
        }

    def compute_energy(self, oth):

        """Computes the energy between two H-bonding groups"""

        if oth.atoms["H"] is not None:
            q_1 = 0.42
            q_2 = 0.2
            dim_f = 332
            dist_on = self.atoms["O"].compute_distance(oth.atoms["N"])
            dist_ch = self.atoms["C"].compute_distance(oth.atoms["H"])
            dist_oh = self.atoms["O"].compute_distance(oth.atoms["H"])
            dist_cn = self.atoms["C"].compute_distance(oth.atoms["N"])
            energy = q_1*q_2*((1/dist_on)+(1/dist_ch)-(1/dist_oh)-(1/dist_cn))*dim_f
        else:
            energy = 0
        return energy

    def add_h(self, oth):

        """Adds hydrogen atoms"""

        v_n = self.atoms["N"].pos_vector()
        v_co = oth.atoms["C"].vector(oth.atoms["O"])
        dist_co = oth.atoms["C"].compute_distance(oth.atoms["O"])
        v_h = vectors.compute_frac_vect(v_co, dist_co)
        v_nh = vectors.compute_sum_vect(v_n, v_h)
        h_atom = ["H", self.name, self.chain, self.number, v_nh]
        self.atoms["H"] = Atom(h_atom, hydrogen=True)

    def compute_tco(self, oth):

        """Computes the tco"""

        v_i = self.atoms["C"].vector(self.atoms["O"])
        v_j = oth.atoms["C"].vector(oth.atoms["O"])
        angle = math.cos(angles.compute_angle(v_i, v_j))
        self.angles["TCO"] = angle

    def compute_kappa(self, oth_b, oth_a):

        """Computes kappa"""

        v_i = self.atoms["CA"].pos_vector()
        v_j = oth_b.atoms["CA"].pos_vector()
        v_k = oth_a.atoms["CA"].pos_vector()
        v_ij = vectors.compute_diff_vect(v_i, v_j)
        v_ki = vectors.compute_diff_vect(v_k, v_i)
        angle = math.degrees(angles.compute_angle(v_ki, v_ij))
        self.angles["KAPPA"] = angle
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
        self.angles["ALPHA"] = angle
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
        self.angles["PHI"] = angle

    def compute_psi(self, oth):

        """Compute psi"""

        v_i = self.atoms["N"].pos_vector()
        v_j = self.atoms["CA"].pos_vector()
        v_k = self.atoms["C"].pos_vector()
        v_l = oth.atoms["N"].pos_vector()
        angle = angles.compute_dihedral_angle(v_i, v_j, v_k, v_l)
        self.angles["PSI"] = angle

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
        "TRP":"W", "TYR":"Y", "UNK":"K"
    }

    def __init__(self, line, pdb=True, hydrogen=False):

        if hydrogen:
            name, aa_name, chain, number, h_coord = line
            self.atom_name = name
            self.aa_name = aa_name
            self.aa_chain = chain
            self.aa_nb = number
            self.x_coord = h_coord[0]
            self.y_coord = h_coord[1]
            self.z_coord = h_coord[2]
        else:
            if pdb:
                self.atom_name = line[12:16].strip()
                if line[17:20] in self.code:
                    self.aa_name = self.code[line[17:20]]
                else:
                    self.aa_name = self.code["UNK"]
                self.aa_chain = line[21:22]
                self.aa_nb = int(line[22:26])
                self.x_coord = float(line[30:38])
                self.y_coord = float(line[38:46])
                self.z_coord = float(line[46:54])
            else:
                line = line.split()
                self.atom_name = line[3]
                if line[5] in self.code:
                    self.aa_name = self.code[line[5]]
                else:
                    self.aa_name = self.code["UNK"]
                self.aa_chain = line[6]
                self.aa_nb = int(line[8])
                self.x_coord = float(line[10])
                self.y_coord = float(line[11])
                self.z_coord = float(line[12])

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

def add_hydrogen(chains, verbose):

    """Adds hydrogen atoms"""

    if verbose:
        print("Adding H")
    for chain in chains:
        indices = [res.number for res in chain]
        nb_res = len(chain)
        for i in range(nb_res):
            if chain[i].name != "P":
                if chain[i].number-1 in indices:
                    chain[i].add_h(chain[i-1])
                else:
                    chain[i].atoms["H"] = None
            else:
                chain[i].atoms["H"] = None
    if verbose:
        print("Done")

def store_metadata(lines, pdb):

    """Store metadata info"""

    protein_info = {
        "header_pdb":"",
        "organism":"",
        "molecule":"",
        "authors":[],
        "ss_bonds":[]
    }
    if pdb:
        for line in lines:
            if line[0:6] == 'HEADER':
                protein_info["header_pdb"] = line[10:].strip()
            elif '2 ORGANISM_SCIENTIFIC: ' in line:
                protein_info["organism"] += line[9:].strip()
            elif '2 MOLECULE: ' in line:
                protein_info["molecule"] += line[9:].strip()
            elif line[0:6] == 'AUTHOR':
                protein_info["authors"].append(line[10:].strip().split(","))
            elif line[0:6] == 'SSBOND':
                ssbond = [[line[15], float(line[16:21])], [line[29], float(line[30:35])]]
                protein_info["ss_bonds"].append(ssbond)
    protein_info["authors"] = ','.join(sum(protein_info["authors"], []))
    return protein_info

def store_atoms(lines, pdb):

    """Store atoms info"""

    residues = []
    chains = []
    res_nb = 0
    current_chain = ""
    atoms = []
    for line in lines:
        atom = Atom(line, pdb)
        if atom.aa_nb != res_nb:
            if len(atoms) > 0:
                res = Residue(atoms)
                if res.chain != current_chain:
                    if len(residues) > 0:
                        chains.append(residues)
                        residues = []
                    current_chain = res.chain
                residues.append(res)
                atoms = []
            res_nb = atom.aa_nb
        atoms.append(atom)
    res = Residue(atoms)
    residues.append(res)
    chains.append(residues)
    return chains

def read_protein_file(lines, pdb):

    """Read lines"""

    metadata = []
    atoms = []
    for line in lines:
        if not line.startswith('ATOM') and not line.startswith('HETATM'):
            metadata.append(line)
        elif line.startswith('ATOM'):
            atoms.append(line)
    try:
        checks.check_empty_protein(atoms)
        protein_info = store_metadata(metadata, pdb)
        chains = store_atoms(atoms, pdb)
        return protein_info, chains
    except AssertionError as error:
        sys.exit(error)

def main():

    """Main function"""

    # Get the arguments
    method_help = (
        "the secondary structure assignment method to use: "
        "ssam, dssp, ssamcompare, or dsspcompare"
    )
    parser = argparse.ArgumentParser()
    parser.add_argument("ssam", type=str, help=method_help)
    parser.add_argument("-v", "--verbose", action="store_true", help="increase output verbosity")
    parser.add_argument("-i", type=str, help="the input file, either a .pdb or .cif")
    parser.add_argument("-o", type=str, help="the output file, a .dssp file")
    parser.add_argument("-hy", "--hydrogen", action="store_true", help="add hydrogen atoms")
    args = parser.parse_args()

    # Read the file
    try:
        checks.check_file_extension(args.i)
        with open(args.i, "r") as input_file:
            if args.verbose:
                print("Reading file")
            lines = input_file.readlines()
            if args.i.lower().endswith(".pdb"):
                protein_info, chains = read_protein_file(lines, pdb=True)
            else:
                protein_info, chains = read_protein_file(lines, pdb=False)
            for chain in chains:
                checks.check_residues_order(chain)
            if args.verbose:
                print("Done")
    except AssertionError as error:
        sys.exit(error)
    except FileNotFoundError:
        sys.exit("{} does not exist".format(args.i))

    # Remove incomplete residues
    for chain in chains:
        res_to_remove = []
        for res in chain:
            try:
                checks.check_residues_completeness(res)
            except AssertionError as error:
                print(error)
                res_to_remove.append(res)
        for res in res_to_remove:
            chain.remove(res)
    chains = [chain for chain in chains if len(chain) > 0]
    
    # Exit if the protein has valid residues
    try:
        checks.check_valid_residues(chains)
    except AssertionError as error:
        sys.exit(error)

    # Add hydrogen atoms if specified
    if args.hydrogen:
        add_hydrogen(chains, args.verbose)
    else:
        try:
            checks.check_hydrogen(chains)
        except AssertionError as error:
            print(error)
            add_hydrogen(chains, args.verbose)

    # SSAM
    try:
        checks.check_method(args.ssam)
        if args.ssam == "dssp":
            # DSSP
            dssp.dssp(args.o, protein_info, chains, args.verbose)
        elif args.ssam == "dsspcompare":
            # DSSP
            dssp.dssp(args.o, protein_info, chains, args.verbose)
            # DSSP compare
            dssp.dssp_compare(args.i, args.o, chains, args.verbose)
        elif args.ssam == "ssam":
            # DSSP
            dssp.dssp(args.o, protein_info, chains, args.verbose)
        elif args.ssam == "ssamcompare":
            # DSSP
            dssp.dssp(args.o, protein_info, chains, args.verbose)
            # DSSP compare
            dssp.dssp_compare(args.i, args.o, chains, args.verbose)

    except AssertionError as error:
        sys.exit(error)

if __name__ == '__main__':
    main()
