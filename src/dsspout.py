"""Output text module.

This module generates the text to write to the dssp file. For instance, it
writes the information about each residue.
"""

from datetime import date

def out_pdb_info(outfile, pdb_info):

    """Write pdb info"""

    today = date.today()
    header_pdb, organism, molecule, authors = pdb_info
    header_dssp = (
        "==== Secondary Structure Definition by the program DSSP, "
        "CMBI version by M.L. Hekkelman/2010-10-21 ==== DATE="
        "{:<18s}.\n".format(today.strftime("%Y-%m-%d"))
    )
    article_ref = (
        "REFERENCE W. KABSCH AND C.SANDER, BIOPOLYMERS 22 (1983) 2577-2637     "
        "                                                         .\n"
    )
    outfile.write(header_dssp)
    outfile.write(article_ref)
    outfile.write("{:<127s}.\n".format(header_pdb))
    outfile.write("{:<127s}.\n".format(molecule))
    outfile.write("{:<127s}.\n".format(organism))
    outfile.write("AUTHOR    {:<117s}.\n".format(authors))

def out_residues(out_file, residues):

    """Write information about each residue"""

    i = 1
    prev_res = residues[0].number-1
    summary_dssp = (
        "  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O"
        "    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA"
    )
    out_file.write(summary_dssp+"\n")
    for res in residues:
        if res.number - prev_res != 1:
            line = (
                "{:>5d}{:>5s} {:>1s} {:>1s}  "
                "{:>1s} {:>1s}{:>1s}{:>1s}"
                "{:>1s}{:>1s}"
                "{:62s}"
                "{:6.3f}{:6.1f}{:6.1f}{:6.1f}{:6.1f}"
                "{:7.1f}{:7.1f}{:7.1f}"
                "\n".format(i, "", "", "!",
                "", "", "", "", 
                "", "",
                "",
                0, 360, 360, 360, 360,
            	0, 0, 0)
            )
            out_file.write(line)
            i += 1
        line = (
            "{:>5d}{:>5d} {:>1s} {:>1s}  "
            "{:>1s} {:>1s}{:>1s}{:>1s}"
            "{:>1s}{:>1s}"
            "{:62s}"
            "{:6.3f}{:6.1f}{:6.1f}{:6.1f}{:6.1f}"
            "{:7.1f}{:7.1f}{:7.1f}"
            "\n".format(i, res.number, res.chain, res.name,
            res.struct["STRC"], res.struct["3"], res.struct["4"], res.struct["5"], 
            res.struct["BEND"], res.struct["CHR"],
            "",
            res.angles["TCO"], res.angles["KAPPA"], res.angles["ALPHA"], 
			res.angles["PHI"], res.angles["PSI"],
			res.atoms["CA"].x_coord, res.atoms["CA"].y_coord, res.atoms["CA"].z_coord)
        )
        out_file.write(line)
        i += 1
        prev_res = res.number
