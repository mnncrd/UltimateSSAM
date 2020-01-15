"""Checking module.

This module performs checks.
"""

def check_nb_res(struct_ssam, struct_comp):

    """Checks if there is the same number of residues in both files"""

    assert len(struct_ssam) == len(struct_comp), "Not the same number of residues in both files"

def check_file_extension(filename):

    """Checks if the file has a.pdb or .cif extension"""

    assert (filename.lower().endswith(".pdb") or filename.lower().endswith(".cif")), (
        "{} is not a .pdb or .cif (PDBx) file".format(filename)
    )

def check_residues_order(residues):

    """Checks if the residues are stored in sequential order"""

    res_nb = [residue.number for residue in residues]
    ord_res_nb = sorted(res_nb)
    assert res_nb == ord_res_nb, "The residues are not in sequential order"

def check_method(ssam_method):

    """Checks if the secondary structure assignment method is available"""

    ssam_methods = ["dssp", "ssam", "dsspcompare", "ssamcompare"]
    assert ssam_method in ssam_methods, (
        "{} is not a secondary structure assignment method available in "
        "UltimateSSAM".format(ssam_method)
    )

def check_residues_completeness(res):

    """Checks if residue is complete"""

    code = {
        "A":"ALA", "C":"CYS", "D":"ASP", "E":"GLU", "F":"PHE", "G":"GLY",
        "H":"HIS", "I":"ILE", "K":"LYS", "L":"LEU", "M":"MET", "N":"ASN",
        "P":"PRO", "Q":"GLN", "R":"ARG", "S":"SER", "T":"THR", "V":"VAL",
        "W":"TRP", "Y":"TYR"
    }
    missing_CA = res.atoms["CA"] is None
    missing_C = res.atoms["C"] is None
    missing_O = res.atoms["O"] is None
    missing_N = res.atoms["N"] is None
    complete_res = not(missing_CA or missing_C or missing_O or missing_N)
    assert complete_res, (
        "Ignoring incomplete residue {:s} ({:d})".format(code[res.name], res.number)
    )

def check_valid_residues(chains):

    """Checks has valid residues"""

    assert len(chains) > 0, "The protein has no valid complete residues"

def check_empty_protein(atoms):

    """Checks if the protein is empty"""

    assert len(atoms) > 0, "The protein is empty"

def check_hydrogen(chains):

    """Checks if hydrogens are present"""

    nb_H = len([res.atoms["H"] for chain in chains for res in chain if res.atoms["H"] is not None])
    assert nb_H > 0, "No hydrogen atoms present, UltimateSSAM will add them"
