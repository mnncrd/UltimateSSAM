"""Hydrogen bonds module.

This module finds hydrogen bonds in the protein.
"""

def assign_hbonds(residues):

    """Assigns h_bonds to residues"""

    nb_res = len(residues)
    for i in range(nb_res):
        all_hb = {}
        for j in range(nb_res):
            if i != j:
                energy = residues[i].compute_energy(residues[j])
                if energy > -3.5:
                    all_hb[residues[j].number - residues[i].number] = energy
        min_hb = min(all_hb, key=all_hb.get)
        residues[i].bonds["O1"] = min_hb
        residues[i].bonds["VO1"] = all_hb[min_hb]
        del all_hb[min_hb]
        min_hb = min(all_hb, key=all_hb.get)
        residues[i].bonds["O2"] = min_hb
        residues[i].bonds["VO2"] = all_hb[min_hb]
    for i in range(1, nb_res):
        all_hb = {}
        for j in range(nb_res):
            if i != j:
                energy = residues[j].compute_energy(residues[i])
                if energy > -3.5:
                    all_hb[residues[j].number - residues[i].number] = energy
        min_hb = min(all_hb, key=all_hb.get)
        residues[i].bonds["N1"] = min_hb
        residues[i].bonds["VN1"] = all_hb[min_hb]
        del all_hb[min_hb]
        min_hb = min(all_hb, key=all_hb.get)
        residues[i].bonds["N2"] = min_hb
        residues[i].bonds["VN2"] = all_hb[min_hb]


def nb_hbonds(para_hbonds, anti_hbonds, residues):

    """Finds the number of hydrogen bonda"""

    nb_para_hbonds = len(para_hbonds)
    nb_anti_hbonds = len(anti_hbonds)
    h_bond_type = []
    for res in residues:
        if res.bonds["VO1"] < -0.5:
            h_bond_type.append(res.bonds["O1"])
        if res.bonds["VO2"] < -0.5:
            h_bond_type.append(res.bonds["O2"])
    nb_h_bond_type = [0]*11
    for i in range(11):
        nb_h_bond_type[i] = h_bond_type.count(i-5)
    nb_tot = len(h_bond_type)
    h_bonds = [nb_tot, nb_para_hbonds, nb_anti_hbonds]
    h_bonds.extend(nb_h_bond_type)
    return h_bonds
