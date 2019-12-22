"""Output text module.

This module generates the text to write to the dssp file. For instance, it
writes the information about each residue.
"""

def out_residues(out_file, residues):

    """Write information about each residue"""

    summary_dssp = (
        "  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O"
        "    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA"
    )
    out_file.write(summary_dssp+"\n")
    for i, res in enumerate(residues):
        line = (
            "{:>5d}{:>5d} {:>1s} {:>1s}  "
            "{:>1s} {:>1s}{:>1s}{:>1s}"
            "{:>1s}{:>1s}"
            "{:62s}"
            "{:6.3f}{:6.1f}{:6.1f}{:6.1f}{:6.1f}"
            "{:7.1f}{:7.1f}{:7.1f}"
            "\n".format(i+1, res.number, res.chain, res.name,
            "", res.struct["3"], res.struct["4"], res.struct["5"], 
            res.struct["BEND"], res.struct["CHR"],
            "",
            res.angles["TCO"], res.angles["KAPPA"], res.angles["ALPHA"], 
			res.angles["PHI"], res.angles["PSI"],
			res.atoms["CA"].x_coord, res.atoms["CA"].y_coord, res.atoms["CA"].z_coord)
        )
        out_file.write(line)
