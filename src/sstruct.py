"""Secondary structures assignment module.

This module checks if a residue meets the requirement for every secondary structure.
For instance, parallel bridges.
"""

def anti_bridge(residues):

    """Finds antiparallel bridges"""

    abridges = []
    indices = [residue.number for residue in residues]
    nb_res = len(residues)
    for i in range(1, nb_res-4):
        if residues[i].number-1 in indices and residues[i].number+1 in indices:
            for j in range(i+3, nb_res-1):
                if residues[j].number-1 in indices and residues[j].number+1 in indices:
                    energy_1 = residues[i].compute_energy(residues[j])
                    energy_2 = residues[j].compute_energy(residues[i])
                    energy_3 = residues[i-1].compute_energy(residues[j+1])
                    energy_4 = residues[j-1].compute_energy(residues[i+1])
                    if energy_1 < -0.5 and energy_2 < -0.5:
                        bridge = (residues[i], residues[j])
                        residues[i].struct["B"] = True
                        residues[j].struct["B"] = True
                        abridges.append(bridge)
                    elif energy_3 < -0.5 and energy_4 < -0.5:
                        bridge = (residues[i], residues[j])
                        residues[i].struct["B"] = True
                        residues[j].struct["B"] = True
                        abridges.append(bridge)
    return abridges

def para_bridge(residues):

    """Finds parallel bridges"""

    pbridges = []
    indices = [residue.number for residue in residues]
    nb_res = len(residues)
    for i in range(1, nb_res-4):
        if residues[i].number-1 in indices and residues[i].number+1 in indices:
            for j in range(i+3, nb_res-1):
                if residues[j].number-1 in indices and residues[j].number+1 in indices:
                    energy_1 = residues[i-1].compute_energy(residues[j])
                    energy_2 = residues[j].compute_energy(residues[i+1])
                    energy_3 = residues[j-1].compute_energy(residues[i])
                    energy_4 = residues[i].compute_energy(residues[j+1])
                    if energy_1 < -0.5 and energy_2 < -0.5:
                        bridge = (residues[i], residues[j])
                        residues[i].struct["B"] = True
                        residues[j].struct["B"] = True
                        pbridges.append(bridge)
                    elif energy_3 < -0.5 and energy_4 < -0.5:
                        bridge = (residues[i], residues[j])
                        residues[i].struct["B"] = True
                        residues[j].struct["B"] = True
                        pbridges.append(bridge)
    return pbridges
