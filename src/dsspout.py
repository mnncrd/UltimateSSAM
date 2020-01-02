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
        "version by M. Curaudeau/2020-01-02 ==== DATE="
    )
    article_ref = ("REFERENCE W. KABSCH AND C.SANDER, BIOPOLYMERS 22 (1983) 2577-2637")
    outfile.write("{:<s}{:<25s}.\n".format(header_dssp, today.strftime("%Y-%m-%d")))
    outfile.write("{:<127s}.\n".format(article_ref))
    outfile.write("{:<127s}.\n".format(header_pdb))
    outfile.write("{:<127s}.\n".format(molecule))
    outfile.write("{:<127s}.\n".format(organism))
    outfile.write("AUTHOR    {:<117s}.\n".format(authors))

def out_histogram(out_file, structures):

    """Write histogram"""

    helices, anti_ladders, para_ladders, sheets = structures
    # Header
    for i in range(1, 31):
        out_file.write("{:>3d}".format(i))
    out_file.write("{:>4s}{:<33s}.\n".format("", " *** HISTOGRAMS OF *** "))
    # Residues per alpha helix
    alpha_helices = helices[1]
    res_alpha_helices = [len(helix) for helix in alpha_helices]
    for i in range(1, 31):
        out_file.write("{:>3d}".format(res_alpha_helices.count(i)))
    out_file.write("{:>4s}{:<33s}.\n".format("", "RESIDUES PER ALPHA HELIX"))
    # Parallel bridges per ladder
    bridges_para_ladders = [len(ladder) for ladder in para_ladders]
    for i in range(1, 31):
        out_file.write("{:>3d}".format(bridges_para_ladders.count(i)))
    out_file.write("{:>4s}{:<33s}.\n".format("", "PARALLEL BRIDGES PER LADDER"))
    # Antiparallel bridges per ladder
    bridges_anti_ladders = [len(ladder) for ladder in anti_ladders]
    for i in range(1, 31):
        out_file.write("{:>3d}".format(bridges_anti_ladders.count(i)))
    out_file.write("{:>4s}{:<33s}.\n".format("", "ANTIPARALLEL BRIDGES PER LADDER"))
    # Ladders per sheet
    ladder_sheets = [len(sheet) for sheet in sheets]
    for i in range(1, 31):
        out_file.write("{:>3d}".format(ladder_sheets.count(i)))
    out_file.write("{:>4s}{:<33s}.\n".format("", "LADDERS PER SHEET"))


def out_residues(out_file, residues):

    """Write information about each residue"""

    i = 1
    current_chain = residues[0].chain
    prev_res = residues[0].number-1
    summary_dssp = (
        "  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O"
        "    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA"
    )
    out_file.write(summary_dssp+"\n")
    for res in residues:
        if res.chain != current_chain:
            line = (
                "{:>5d}{:>8s}{:>2s}{:>10s}"
                "{:4d}{:4d}{:1s}{:51s}"
                "{:6.3f}{:6.1f}{:6.1f}{:6.1f}{:6.1f}"
                "{:7.1f}{:7.1f}{:7.1f}"
                "\n".format(i, "", "!*", "",
                            0, 0, "", "",
                            0, 360, 360, 360, 360,
                            0, 0, 0
                            )
            )
            out_file.write(line)
            i += 1
            current_chain = res.chain
        elif res.number - prev_res != 1:
            line = (
                "{:>5d}{:>8s}{:>1s}{:>11s}"
                "{:4d}{:4d}{:1s}{:51s}"
                "{:6.3f}{:6.1f}{:6.1f}{:6.1f}{:6.1f}"
                "{:7.1f}{:7.1f}{:7.1f}"
                "\n".format(i, "", "!", "",
                            0, 0, "", "",
                            0, 360, 360, 360, 360,
                            0, 0, 0
                            )
            )
            out_file.write(line)
            i += 1
        line = (
            "{:>5d}{:>5d} {:>1s} {:>1s}  "
            "{:>1s} {:>1s}{:>1s}{:>1s}"
            "{:>1s}{:>1s}"
            "{:1s}{:1s}{:4d}{:4d}{:1s}"
            "{:51s}"
            "{:6.3f}{:6.1f}{:6.1f}{:6.1f}{:6.1f}"
            "{:7.1f}{:7.1f}{:7.1f}"
            "\n".format(i, res.number, res.chain, res.name,
                        res.struct["STRC"], res.struct["3"], res.struct["4"], res.struct["5"],
                        res.struct["BEND"], res.struct["CHR"],
                        res.struct["LAD1"], res.struct["LAD2"],
                        res.struct["BP1"], res.struct["BP2"],
                        res.struct["SHEET"],
                        "",
                        res.angles["TCO"], res.angles["KAPPA"], res.angles["ALPHA"],
                        res.angles["PHI"], res.angles["PSI"],
                        res.atoms["CA"].x_coord, res.atoms["CA"].y_coord, res.atoms["CA"].z_coord
                        )
        )
        out_file.write(line)
        i += 1
        prev_res = res.number
