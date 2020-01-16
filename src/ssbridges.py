"""SS-bridges module.

This module deals with SS-bridges/SS-bonds. For instance, it counts the number
of intra and inter-SS-bonds. It also assigns SS-bridges to residues.
"""

def count_ss_bonds(ss_bonds):

    """Counts the number of intra and inter-SS-bonds.

    For each bond in the given list, checks whether the two atoms are in the
    same chain.

    Args:
        ss_bonds: A list of bonds.

    Returns:
        A tuple contianing the number of intra-SS-bonds.
    """

    intra = 0
    inter = 0
    for bond in ss_bonds:
        if bond[0][0] == bond[1][0]:
            intra += 1
        else:
            inter += 1
    nb_ss = (intra, inter)
    return nb_ss

def assign_ss_bonds(residues, ss_bonds):

    """Assigns SS-bond to each residue in an SS-bridge.

    For each Residue instance in the given residues list, change the residue
    name to reflect whether said residue is involved in an SS-bridge.

    Args:
        residues: A list of Residue instances.
        ss_bonds: A list of bonds.
    """

    alphabet = {i:chr(65+i).lower() for i in range(26)}
    for res in residues:
        for i, bond in enumerate(ss_bonds):
            if res.chain == bond[0][0] and res.number == bond[0][1]:
                res.name = alphabet[i%26]
            elif res.chain == bond[1][0] and res.number == bond[1][1]:
                res.name = alphabet[i%26]
