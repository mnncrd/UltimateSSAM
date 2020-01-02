"""SS-bridges module.

This module deals with SS-bridges. For instance, it assign ss-bridges.
"""

def count_ss_bonds(ss_bonds):

    intra = 0
    inter = 0
    for bond in ss_bonds:
        if bond[0][0] == bond[1][0]:
            intra += 1
        else:
            inter += 1
    nb_ss = (intra, inter)
    return(nb_ss)

def assign_ss_bonds(residues, ss_bonds):

    """Assign SS-bonds"""

    alphabet = {i:chr(65+i).lower() for i in range(26)}
    for res in residues:
        for i, bond in enumerate(ss_bonds):
            if res.chain == bond[0][0] and res.number == bond[0][1]:
                res.name = alphabet[i%26]
            elif res.chain == bond[1][0] and res.number == bond[1][1]:
                res.name = alphabet[i%26]