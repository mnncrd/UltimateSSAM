"""DSSP module.

This module does DSSP.
"""

# Import libraries
import os
import sys

# Import modules
import dsspout
import sstruct
import ssbridges
import hbonds
import checks

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

def dssp(outfile, protein_info, chains, verbose):

    """Does DSSP"""

	# structures = [Helices, Anti-parallel ladders, Parallel ladders, Sheets]
    structures = [[[], [], []], [], [], []]
    h_bonds = [0]*14
    for chain in chains:
        if verbose:
            print("Chain {:s}".format(chain[0].chain))
            print("Computing angles")
        find_angles(chain)
        if verbose:
            print("Done")
            print("Assign hydrogen bonds")
        hbonds.assign_hbonds(chain)
        if verbose:
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
        if verbose:
            print("Done")
    if verbose:
        print("Write .dssp file")
    write_dssp_file(outfile, protein_info, chains, structures, h_bonds)
    if verbose:
        print("Done")

def stats_dssp(struct_ssam, struct_comp):

    """Performs stats"""

    stats = {}
    spl_strct = {"G":"H", "H":"H", "I":"H", "B":"E", "E":"E", "S":"C", "T":"C", " ":"C"}
    nb_match = 0
    nb_res = len(struct_ssam)
    for s_ssam, s_comp in zip(struct_ssam, struct_comp):
        if s_ssam == s_comp:
            nb_match += 1
    stats["acc"] = nb_match/nb_res
    nb_match_3 = {"H":0, "E":0, "C":0}
    for s_ssam, s_comp in zip(struct_ssam, struct_comp):
        if spl_strct[s_ssam] == spl_strct[s_comp]:
            nb_match_3[spl_strct[s_ssam]] += 1
    nb_h = len([x for x in struct_comp if spl_strct[x] == "H"])
    nb_e = len([x for x in struct_comp if spl_strct[x] == "E"])
    nb_c = len([x for x in struct_comp if spl_strct[x] == "C"])
    if nb_h > 0:
        stats["acch"] = nb_match_3["H"]/nb_h
    else:
        stats["acch"] = 1.0
    if nb_e > 0:
        stats["acce"] = nb_match_3["E"]/nb_e
    else:
        stats["acce"] = 1.0
    if nb_c > 0:
        stats["accc"] = nb_match_3["C"]/nb_c
    else:
        stats["accc"] = 1.0
    return stats

def dssp_compare(infile, outfile, outcompare, chains, verbose):

    """Launches a comparison to DSSP"""

    dssp_file = "{:s}_dssp.dssp".format(outfile[:-5])
    cmd = "mkdssp -i {:s} -o {:s}".format(infile, dssp_file)
    if verbose:
        print("Launch a comparison to DSSP")
    os.system(cmd)
    struct_ssam = [res.struct["STRC"] for chain in chains for res in chain]
    struct_comp = []
    with open(dssp_file, "r") as input_file:
        while not '#' in next(input_file):
            pass
        for line in input_file:
            if '!' not in line:
                struct_comp.append(line[16])
    try:
        checks.check_nb_res(struct_ssam, struct_comp)
        stats = stats_dssp(struct_ssam, struct_comp)
        with open(outcompare, "a") as output_file:
            line = (
                "{:s},{:4f},{:4f},{:4f},{:4f}"
                "\n".format(infile, stats["acc"], stats["acch"], stats["acce"], stats["accc"])
            )
            output_file.write(line)
        if verbose:
            print("Done")
    except AssertionError as error:
        sys.exit(error)
