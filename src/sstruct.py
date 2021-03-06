"""Secondary structures assignment module.

This module checks if a residue meets the requirement for every secondary structure.
For instance, parallel bridges.
"""

def structure_to_print(residues):

    """Finds which structure to print"""

    for res in residues:
        if res.struct["S"]:
            res.struct["STRC"] = "S"
        if res.struct["T"]:
            res.struct["STRC"] = "T"
        if res.struct["I"]:
            res.struct["STRC"] = "I"
        if res.struct["G"]:
            res.struct["STRC"] = "G"
        if res.struct["B"]:
            res.struct["STRC"] = "B"
        if res.struct["E"]:
            res.struct["STRC"] = "E"
        if res.struct["H"]:
            res.struct["STRC"] = "H"

def assign_sheets(sheets):

    """Assigns sheet name"""

    alphabet = {i:chr(65+i) for i in range(26)}
    for i, sht in enumerate(sheets):
        for lad in sht:
            for res in lad:
                res.struct["SHEET"] = alphabet[i%26]

def sheet(para_ladders, anti_ladders):

    """Finds sheets"""

    ladders = para_ladders + anti_ladders
    res_in_ladders = []
    for ladder in ladders:
        res = []
        if len(ladder) > 0:
            for bridge in ladder:
                res.extend([bridge[0], bridge[1]])
            res_in_ladders.append(res)
    res_in_sheets = list(res_in_ladders)
    sheets = []
    l_count = len(res_in_sheets)
    if l_count > 1:
        l_count = len(res_in_sheets)
        for i in range(l_count-1):
            for j in range(i+1, l_count):
                if len(set(res_in_sheets[i]) & set(res_in_sheets[j])) > 0:
                    sht = [res_in_sheets[i], res_in_sheets[j]]
                    sheets.append(sht)
                    res_in_sheets[i] = list(set(res_in_sheets[i]) | set(res_in_sheets[j]))
                    res_in_sheets[j] = []
    sheets.extend([[lad] for lad in res_in_ladders if lad in res_in_sheets])
    res_in_sheets = [sorted(lad, key=lambda res: res.number)
                     for lad in res_in_sheets if len(lad) > 0]
    s_count = len(sheets)
    if s_count > 1:
        for i in range(s_count-1):
            for j in range(i+1, s_count):
                if len(set(sheets[i][-1]) & set(sheets[j][-1])) > 0:
                    sheets[i].append(sheets[j][-1])
                    sheets[j] = [[]]
    sheets = [[sorted(lad, key=lambda res: res.number) for lad in sht if len(sht[0]) > 0]
              for sht in sheets]
    sheets = [sht for sht in sheets if len(sht) > 0]
    assign_sheets(sheets)
    return sheets

def assign_ladder(alphabet, ladders):

    """Assigns ladder name"""

    for i, lad in enumerate(ladders):
        res_in_lad_1, res_in_lad_2 = 0, 0
        for bridge in lad:
            if bridge[0].struct["LAD1"] != "":
                res_in_lad_1 += 1
            if bridge[1].struct["LAD1"] != "":
                res_in_lad_2 += 1
        for bridge in lad:
            if res_in_lad_1 == 0:
                bridge[0].struct["LAD1"] = alphabet[i%26]
                bridge[0].struct["BP1"] = bridge[1].number
            else:
                bridge[0].struct["LAD2"] = alphabet[i%26]
                bridge[0].struct["BP2"] = bridge[1].number
            if res_in_lad_2 == 0:
                bridge[1].struct["LAD1"] = alphabet[i%26]
                bridge[1].struct["BP1"] = bridge[0].number
            else:
                bridge[1].struct["LAD2"] = alphabet[i%26]
                bridge[1].struct["BP2"] = bridge[0].number

def anti_ladder(bridges):

    """Finds ladders"""

    indices_bridges = [(bridge[0].number, bridge[1].number) for bridge in bridges]
    bridge_in_ladder = []
    ladders = []
    for i in indices_bridges:
        idx_first_bridge = indices_bridges.index(i)
        if (i[0]+1, i[1]-1) in indices_bridges:
            idx_second_bridge = indices_bridges.index((i[0]+1, i[1]-1))
            if bridges[idx_first_bridge] not in bridge_in_ladder:
                ladder = [bridges[idx_first_bridge], bridges[idx_second_bridge]]
                bridge_in_ladder.append(bridges[idx_first_bridge])
                bridge_in_ladder.append(bridges[idx_second_bridge])
                bridges[idx_first_bridge][0].struct["E"] = True
                bridges[idx_first_bridge][1].struct["E"] = True
                bridges[idx_second_bridge][0].struct["E"] = True
                bridges[idx_second_bridge][1].struct["E"] = True
                ladders.append(ladder)
            else:
                bridge_in_ladder.append(bridges[idx_second_bridge])
                for ladder in ladders:
                    if bridges[idx_first_bridge] in ladder:
                        bridges[idx_second_bridge][0].struct["E"] = True
                        bridges[idx_second_bridge][1].struct["E"] = True
                        ladder.append(bridges[idx_second_bridge])
        else:
            if bridges[idx_first_bridge] not in bridge_in_ladder:
                ladders.append([bridges[idx_first_bridge]])
    alphabet = {i:chr(65+i) for i in range(26)}
    assign_ladder(alphabet, ladders)
    return ladders

def para_ladder(bridges):

    """Finds ladders"""

    indices_bridges = [(bridge[0].number, bridge[1].number) for bridge in bridges]
    bridge_in_ladder = []
    ladders = []
    for i in indices_bridges:
        idx_first_bridge = indices_bridges.index(i)
        if (i[0]+1, i[1]+1) in indices_bridges:
            idx_second_bridge = indices_bridges.index((i[0]+1, i[1]+1))
            if bridges[idx_first_bridge] not in bridge_in_ladder:
                ladder = [bridges[idx_first_bridge], bridges[idx_second_bridge]]
                bridge_in_ladder.append(bridges[idx_first_bridge])
                bridge_in_ladder.append(bridges[idx_second_bridge])
                bridges[idx_first_bridge][0].struct["E"] = True
                bridges[idx_first_bridge][1].struct["E"] = True
                bridges[idx_second_bridge][0].struct["E"] = True
                bridges[idx_second_bridge][1].struct["E"] = True
                ladders.append(ladder)
            else:
                bridge_in_ladder.append(bridges[idx_second_bridge])
                for ladder in ladders:
                    if bridges[idx_first_bridge] in ladder:
                        bridges[idx_second_bridge][0].struct["E"] = True
                        bridges[idx_second_bridge][1].struct["E"] = True
                        ladder.append(bridges[idx_second_bridge])
        else:
            if bridges[idx_first_bridge] not in bridge_in_ladder:
                ladders.append([bridges[idx_first_bridge]])
    alphabet = {i:chr(65+i).lower() for i in range(26)}
    assign_ladder(alphabet, ladders)
    return ladders

def final_bridges(anti_bridges, para_bridges):

    """Assign a bridge"""

    bridges = anti_bridges + para_bridges
    for bridge in bridges:
        bridge[0].struct["B"] = True
        bridge[1].struct["B"] = True
    return anti_bridges, para_bridges

def not_in_helix(bridges, helices):

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
    ahbonds = []
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
                        hbond1 = (residues[i].number, residues[j].number)
                        hbond2 = (residues[j].number, residues[i].number)
                        bridge = (residues[i], residues[j])
                        residues[i].struct["B"] = True
                        residues[j].struct["B"] = True
                        ahbonds.append(hbond1)
                        ahbonds.append(hbond2)
                        abridges.append(bridge)
                    elif energy_3 < -0.5 and energy_4 < -0.5:
                        hbond1 = (residues[i-1].number, residues[j+1].number)
                        hbond2 = (residues[j-1].number, residues[i+1].number)
                        bridge = (residues[i], residues[j])
                        residues[i].struct["B"] = True
                        residues[j].struct["B"] = True
                        ahbonds.append(hbond1)
                        ahbonds.append(hbond2)
                        abridges.append(bridge)
    return abridges, list(set(ahbonds))

def para_bridge(residues):

    """Finds parallel bridges"""

    pbridges = []
    phbonds = []
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
                        hbond1 = (residues[i-1].number, residues[j].number)
                        hbond2 = (residues[j].number, residues[i+1].number)
                        bridge = (residues[i], residues[j])
                        residues[i].struct["B"] = True
                        residues[j].struct["B"] = True
                        phbonds.append(hbond1)
                        phbonds.append(hbond2)
                        pbridges.append(bridge)
                    elif energy_3 < -0.5 and energy_4 < -0.5:
                        hbond1 = (residues[j-1].number, residues[i].number)
                        hbond2 = (residues[i].number, residues[j+1].number)
                        bridge = (residues[i], residues[j])
                        residues[i].struct["B"] = True
                        residues[j].struct["B"] = True
                        phbonds.append(hbond1)
                        phbonds.append(hbond2)
                        pbridges.append(bridge)
    return pbridges, list(set(phbonds))

def sort_helices(helices, min_len):

    """Sort residues by number and filter helices smaller than min_len"""

    helices = [sorted(hlx, key=lambda res: res.number) for hlx in helices if len(hlx) >= min_len]
    return helices

def not_in_bridge(helices, bridges):

    """Remove residues that are in a bridge from helices"""

    h_count = len(helices)
    b_count = len(bridges)
    if h_count > 0 and b_count > 0:
        for i in range(h_count):
            for j in range(b_count):
                if len(set(helices[i]) & set(bridges[j])) > 0:
                    helices[i] = set(helices[i]) - set(bridges[j])
    return helices

def assign_helix(helices, n_val):

    """Assigns a helix to residues in helices and a turn to residues in too
    short helices"""

    h_d = {3:"G", 4:"H", 5:"I"}
    for helice in helices:
        if len(helice) >= n_val:
            for res in helice:
                res.struct[h_d[n_val]] = True
        else:
            for res in helice:
                res.struct["T"] = True

def which_helix(helices_1, helices_2):

    """Removes residues that are in two helices from the lower priority helix"""

    h1_count = len(helices_1)
    h2_count = len(helices_2)
    if h1_count > 0 and h2_count > 0:
        for i in range(h1_count):
            for j in range(h2_count):
                if len(set(helices_1[i]) & set(helices_2[j])) > 0:
                    helices_1[i] = set(helices_1[i]) - set(helices_2[j])
    return helices_1

def final_helices(g_helices, h_helices, i_helices, anti_bridges, para_bridges):

    """Return helices"""

    i_helices = which_helix(i_helices, g_helices)
    i_helices = which_helix(i_helices, h_helices)
    g_helices = which_helix(g_helices, h_helices)
    i_helices = not_in_bridge(i_helices, anti_bridges)
    i_helices = not_in_bridge(i_helices, para_bridges)
    g_helices = not_in_bridge(g_helices, anti_bridges)
    g_helices = not_in_bridge(g_helices, para_bridges)
    assign_helix(g_helices, 3)
    assign_helix(h_helices, 4)
    assign_helix(i_helices, 5)
    g_helices = sort_helices(g_helices, 3)
    h_helices = sort_helices(h_helices, 4)
    i_helices = sort_helices(i_helices, 5)
    helices = [g_helices, h_helices, i_helices]
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
    if h_count > 1:
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
                if residues[i].struct[str(n_val)] != "<":
                    residues[i].struct[str(n_val)] = ">"
                else:
                    residues[i].struct[str(n_val)] = "X"
                if residues[i+n_val].struct[str(n_val)] != ">":
                    residues[i+n_val].struct[str(n_val)] = "<"
                else:
                    residues[i+n_val].struct[str(n_val)] = "X"
                for j in range(1, n_val):
                    residues[i+j].struct["T"] = True
                    if residues[i+j].struct[str(n_val)] == " ":
                        residues[i+j].struct[str(n_val)] = str(n_val)
    return nturns
