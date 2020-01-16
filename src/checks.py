"""Checking module.

This module performs checks. For instance, it checks if the the file extension
is valid, if the residues are stored in sequential order, or if a residue is
complete.
"""

def check_nb_res(struct_ssam, struct_comp):

    """Checks the number of residues in both files.

    Checks if there is the same number of residues in both files.

    Args:
        struct_ssam: A list of Residue instances.
        struct_comp: A list of Residue instances.

    Raises:
        AssertionError: When the number of Residue instances is different
        between the two lists.
    """

    assert len(struct_ssam) == len(struct_comp), "Not the same number of residues in both files"

def check_file_extension(filename):

    """Checks the file extension.

    Checks if the file has a .pdb or .cif extension.

    Args:
        filename: A string.

    Raises:
        AssertionError: When the file name does not end with .pdb or .cif.
    """

    assert (filename.lower().endswith(".pdb") or filename.lower().endswith(".cif")), (
        "{} is not a .pdb or .cif (PDBx) file".format(filename)
    )

def check_residues_order(residues):

    """Checks residues order.

    Checks if the residues are stored in sequential order.

    Args:
        residues: A list of Residue instances.

    Raises:
        AssertionError: When the Residue instances are not sorted in ascending
        numbers.
    """

    res_nb = [residue.number for residue in residues]
    ord_res_nb = sorted(res_nb)
    assert res_nb == ord_res_nb, "The residues are not in sequential order"

def check_method(ssam_method):

    """Checks the secondary structure assignment method.

    Checks if the secondary structure assignment method is available in
    UltimateSSAM.

    Args:
        ssam_method: ssam_method.

    Raises:
        AssertionError: When the method is not available in UltimateSSAM.
    """

    ssam_methods = ["dssp", "ssam", "dsspcompare", "ssamcompare"]
    assert ssam_method in ssam_methods, (
        "{} is not a secondary structure assignment method available in "
        "UltimateSSAM".format(ssam_method)
    )

def check_residues_completeness(res):

    """Checks if a residue is complete.

    Checks if there are CA, C, O, and N atoms in the residue.

    Args:
        res: A Residue instance.

    Raises:
        AssertionError: When the Residue instance is incomplete.
    """

    code = {
        "A":"ALA", "C":"CYS", "D":"ASP", "E":"GLU", "F":"PHE", "G":"GLY",
        "H":"HIS", "I":"ILE", "K":"LYS", "L":"LEU", "M":"MET", "N":"ASN",
        "P":"PRO", "Q":"GLN", "R":"ARG", "S":"SER", "T":"THR", "V":"VAL",
        "W":"TRP", "Y":"TYR", "X":"UNK"
    }
    missing_ca = res.atoms["CA"] is None
    missing_c = res.atoms["C"] is None
    missing_o = res.atoms["O"] is None
    missing_n = res.atoms["N"] is None
    complete_res = not(missing_ca or missing_c or missing_o or missing_n)
    assert complete_res, (
        "Ignoring incomplete residue {:s} ({:d})".format(code[res.name], res.number)
    )

def check_valid_residues(chains):

    """Checks if the protein has valid residues.

    Checks if the given list is not empty.

    Args:
        chains: A list of list of Residue instances.

    Raises:
        AssertionError: When the list is empty.
    """

    assert len(chains) > 0, "The protein has no valid complete residues"

def check_empty_protein(atoms):

    """Checks if the protein is empty, i.e. when there is no ATOM lines.

    Checks if the given list is not empty.

    Args:
        atoms: A list of list of Residue instances.

    Raises:
        AssertionError: When the list is empty.
    """

    assert len(atoms) > 0, "The protein is empty"

def check_hydrogen(chains):

    """Checks if hydrogen atoms are present.

    Chhecks the number of Residue instances whith hydrogen atoms present.

    Args:
        chains: A list of list of Residue instances.

    Raises:
        AssertionError: When there is no hydrogen atoms present.
    """

    nb_h = len([res.atoms["H"] for chain in chains for res in chain if res.atoms["H"] is not None])
    assert nb_h > 0, "No hydrogen atoms present, UltimateSSAM will add them"
