"""Secondary structures assignment module.

This module checks if a residue meets the requirement for every secondary structure.
For instance, parallel bridges.
"""

def final_bridges(bridges, helices):

    """Removes bridges that are in an helix"""

    h_count = len(helices)
    b_count = len(bridges)
    if h_count > 0 and b_count > 0:
        for i in range(h_count):
            for j in range(b_count):
                if len(set(helices[i]) & set(bridges[j])) > 0:
                    bridges[j] = ()
    bridges = [bridge for bridge in bridges if len(bridge) == 2]
    return bridges

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

def sort_helices(helices, min_len):

    """Sort residues by number and filter helices smaller than min_len"""

    helices = [sorted(hlx, key=lambda res: res.number) for hlx in helices if len(hlx) >= min_len]
    return helices

def helix(residues, nturns, n_val):

    """Finds helices"""

    helices = []
    indices = [turn.number for turn in nturns]
    nb_res = len(residues)
    for i in range(1, nb_res-n_val):
        hlx = []
        if residues[i].number in indices and residues[i-1].number in indices:
            for j in range(n_val):
                hlx.append(residues[i+j])
            helices.append(hlx)
    h_count = len(helices)
    if h_count > 2:
        for i in range(h_count-1):
            if len(set(helices[i]) & set(helices[i+1])) > 0:
                helices[i+1] = list(set(helices[i]) | set(helices[i+1]))
                helices[i] = []
    helices = sort_helices(helices, n_val)
    return helices

def n_turn(residues, n_val):

    """Finds n-turns"""

    nturns = []
    indices = [residue.number for residue in residues]
    nb_res = len(residues)
    for i in range(nb_res-n_val):
        n_res = 0
        for j in range(1, n_val+1):
            if residues[i].number+j in indices:
                n_res += 1
        if n_res == n_val:
            energy = residues[i].compute_energy(residues[i+n_val])
            if energy < -0.5:
                nturns.append(residues[i])
                residues[i].struct["T"] = True
                if residues[i].struct[str(n_val)] != "<":
                    residues[i].struct[str(n_val)] = ">"
                else:
                    residues[i].struct[str(n_val)] = "X"
                if residues[i+n_val].struct[str(n_val)] != ">":
                    residues[i+n_val].struct[str(n_val)] = "<"
                else:
                    residues[i+n_val].struct[str(n_val)] = "X"
                for j in range(1, n_val):
                    if residues[i+j].struct[str(n_val)] == " ":
                        residues[i+j].struct[str(n_val)] = str(n_val)
    return nturns
