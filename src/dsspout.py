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
