"""SS-bridges module.

This module deals with SS-bridges. For instance, it assign ss-bridges.
"""

def assign_ss_bonds(residues, ss_bonds):

    """Assign SS-bonds"""

    alphabet = {i:chr(65+i).lower() for i in range(26)}
    for res in residues:
        for i, bond in enumerate(ss_bonds):
            if res.chain == bond[0][0] and res.number == bond[0][1]:
                res.name = alphabet[i%26]
            elif res.chain == bond[1][0] and res.number == bond[1][1]:
                res.name = alphabet[i%26]