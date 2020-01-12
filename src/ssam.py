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
import dsspout
import sstruct
import ssbridges
import hbonds

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
        self.number = self.atoms["CA"].aa_nb
        self.name = self.atoms["CA"].aa_name
        self.chain = self.atoms["CA"].aa_chain
        self.angles = {
            "TCO":0.000,
            "KAPPA":360.0,
            "ALPHA":360.0,
            "PHI":360.0,
            "PSI":360.0
        }
        self.struct = {
            "STRC":"",
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
        "TRP":"W", "TYR":"Y"
    }

    def __init__(self, line, pdb, hydrogen=False):

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
                self.aa_name = self.code[line[17:20]]
                self.aa_chain = line[21:22]
                self.aa_nb = int(line[22:26])
                self.x_coord = float(line[30:38])
                self.y_coord = float(line[38:46])
                self.z_coord = float(line[46:54])
            else:
                line = line.split()
                self.atom_name = line[3]
                self.aa_name = self.code[line[5]]
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

def secondary_struct(residues):

    """Assigns secondary structures and find h-bonds in secondary structures"""

    three_turns = sstruct.n_turn(residues, 3)
    four_turns = sstruct.n_turn(residues, 4)
    five_turns = sstruct.n_turn(residues, 5)
    g_helices = sstruct.helix(residues, three_turns, 3)
    h_helices = sstruct.helix(residues, four_turns, 4)
    i_helices = sstruct.helix(residues, five_turns, 5)
    para_bridges, para_hbonds = sstruct.para_bridge(residues)
    anti_bridges, anti_hbonds = sstruct.anti_bridge(residues)
    h_bonds = hbonds.nb_hbonds(para_hbonds, anti_hbonds, residues)
    para_bridges = sstruct.not_in_helix(para_bridges, h_helices)
    anti_bridges = sstruct.not_in_helix(anti_bridges, h_helices)
    anti_bridges, para_bridges = sstruct.final_bridges(anti_bridges, para_bridges)
    para_ladders = sstruct.para_ladder(para_bridges)
    anti_ladders = sstruct.anti_ladder(anti_bridges)
    sheets = sstruct.sheet(para_ladders, anti_ladders)
    helices = sstruct.final_helices(g_helices, h_helices, i_helices, anti_bridges, para_bridges)
    sstruct.structure_to_print(residues)
    structures = (helices, anti_ladders, para_ladders, sheets)
    return structures, h_bonds

def write_dssp_file(filename, protein_info, chains, structures, h_bonds):

    """Create the output .dssp file"""

    file_dssp = open(filename, "w+")
    residues = [res for chain in chains for res in chain]
    nb_chains = len(chains)
    nb_res = len(residues)
    nb_ss = ssbridges.count_ss_bonds(protein_info["ss_bonds"])
    # PDB info
    dsspout.out_protein_info(file_dssp, protein_info)
    # Protein stats
    dsspout.out_stats(file_dssp, nb_res, nb_chains, nb_ss)
    # Accessible surface
    dsspout.out_surface(file_dssp, 0.0)
    # Hydrogen bonds
    dsspout.out_hbonds(file_dssp, h_bonds, nb_res)
    # Histogram
    dsspout.out_histogram(file_dssp, structures)
    # Residues
    ssbridges.assign_ss_bonds(residues, protein_info["ss_bonds"])
    dsspout.out_residues(file_dssp, residues)
    file_dssp.close()

def find_angles(residues):

    """Computes all angles"""

    nb_res = len(residues)
    indices = [residue.number for residue in residues]
    for i in range(nb_res):
        if residues[i].number-1 in indices:
            residues[i].compute_tco(residues[i-1])
            residues[i].compute_phi(residues[i-1])
        if residues[i].number-2 in indices and residues[i].number+2 in indices:
            residues[i].compute_kappa(residues[i-2], residues[i+2])
        if (
                residues[i].number-1 in indices and
                residues[i].number+1 in indices and
                residues[i].number+2 in indices
        ):
            residues[i].compute_alpha(residues[i-1], residues[i+1], residues[i+2])
        if residues[i].number+1 in indices:
            residues[i].compute_psi(residues[i+1])

def read_protein_file(lines, pdb):

    """Stores info"""

    protein_info = {"header_pdb":"", "organism":"", "molecule":"", "authors":[], "ss_bonds":[]}
    protein_info["header_pdb"] = lines[0].strip()
    residues = []
    chains = []
    res_nb = 0
    current_chain = ""
    atoms = []
    for line in lines:
        if '2 ORGANISM_SCIENTIFIC: ' in line:
            protein_info["organism"] += line.strip()
        elif '2 MOLECULE: ' in line:
            protein_info["molecule"] += line.strip()
        elif line[0:6] == 'AUTHOR':
            protein_info["authors"].append(line[10:].strip().split(","))
        elif line[0:6] == 'SSBOND':
            protein_info["ss_bonds"].append([[line[15], float(line[16:21])], [line[29], float(line[30:35])]])
        elif line[0:4] == 'ATOM':
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
    protein_info["authors"] = ','.join(sum(protein_info["authors"], []))
    return protein_info, chains

def check_file_extension(filename):

    """Checks if the file has a.pdb or .cif extension"""

    assert (filename.lower().endswith(".pdb") or filename.lower().endswith(".cif")), (
        "{} is not a .pdb or .cif (PDBx) file".format(filename)
    )

def check_residues_order(residues):

    """Checks if the residues are stored in sequential order"""

    res_nb = [residue.number for residue in residues]
    ord_res_nb = sorted(res_nb)
    assert res_nb == ord_res_nb, "The residues are not in sequential order"

def main():

    """Main function"""

    # Get the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("ssam", help="The secondary structure assignment method to use", type=str)
    parser.add_argument("-v", "--verbose", action="store_true", help="increase output verbosity")
    parser.add_argument("-i", help="The file to assign a secondary structure to", type=str)
    parser.add_argument("-o", type=str)
    parser.add_argument("-hy", "--hydrogen", action="store_true", help="add hydrogen atoms")
    args = parser.parse_args()

    # Read the file
    try:
        check_file_extension(args.i)
        with open(args.i, "r") as input_file:
            if args.verbose:
                print("Reading file")
            lines = input_file.readlines()
            if args.i.lower().endswith(".pdb"):
                protein_info, chains = read_protein_file(lines, pdb=True)
            else:
                protein_info, chains = read_protein_file(lines, pdb=False)
            for chain in chains:
                check_residues_order(chain)
            if args.verbose:
                print("Done")
    except AssertionError as error:
        sys.exit(error)
    except FileNotFoundError:
        sys.exit("{} does not exist".format(args.i))

    # Add hydrogen atoms if specified
    if args.hydrogen:
        if args.verbose:
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

    # DSSP
    # structures = [Helices, Anti-parallel ladders, Parallel ladders, Sheets]
    structures = [[[], [], []], [], [], []]
    h_bonds = [0]*14
    for chain in chains:
        if args.verbose:
            print("Chain {:s}".format(chain[0].chain))
            print("Computing angles")
        find_angles(chain)
        if args.verbose:
            print("Done")
            print("Assign hydrogen bonds")
        hbonds.assign_hbonds(chain)
        if args.verbose:
            print("Done")
            print("Assign secondary structures")
        cur_sec_struct, cur_h_bonds = secondary_struct(chain)
        # Helices
        structures[0] = [x + y for x, y in zip(structures[0], cur_sec_struct[0])]
        # Anti-parallel ladders
        structures[1].extend(cur_sec_struct[1])
        # Parallel ladders
        structures[2].extend(cur_sec_struct[2])
        # Sheets
        structures[3].extend(cur_sec_struct[3])
        # Hydrogen bonds
        h_bonds = [x + y for x, y in zip(h_bonds, cur_h_bonds)]
        if args.verbose:
            print("Done")
    if args.verbose:
        print("Write .dssp file")
    write_dssp_file(args.o, protein_info, chains, structures, h_bonds)
    if args.verbose:
        print("Done")

if __name__ == '__main__':
    main()
