"""Hydrogen bonds module.

This module finds hydrogen bonds in the protein.
"""

def assign_hbonds(residues):

    """Assign h_bonds to residues"""

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
    for i in range(nb_res):
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
