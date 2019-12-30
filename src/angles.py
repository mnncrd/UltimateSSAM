"""Calculation of angles module.

This module computes angles between vectors.
"""

import math
import vectors

def compute_angle(v_i, v_j):

    """Computes the angle between two vectors"""

    scal_prod = vectors.scal_product(v_i, v_j)
    v_i_norm = vectors.vector_norm(v_i)
    v_j_norm = vectors.vector_norm(v_j)
    angle = math.acos(scal_prod/(v_i_norm*v_j_norm))
    return angle

def compute_dihedral_angle(v_i, v_j, v_k, v_l):

    """Computes the dihedral angle between four vectors"""

    v_ji = vectors.compute_diff_vect(v_i=v_j, v_j=v_i)
    v_kj = vectors.compute_diff_vect(v_i=v_k, v_j=v_j)
    v_lk = vectors.compute_diff_vect(v_i=v_l, v_j=v_k)
    cp_ijk = vectors.cross_product(v_ji, v_kj)
    cp_jkl = vectors.cross_product(v_kj, v_lk)
    cp_ijkl = vectors.cross_product(cp_ijk, cp_jkl)
    norm_v_kj = vectors.vector_norm(v_kj)
    unit_v_kj = vectors.compute_frac_vect(v_kj, norm_v_kj)
    angle = math.atan2(vectors.scal_product(cp_ijkl, unit_v_kj),
                       vectors.scal_product(cp_ijk, cp_jkl))
    return math.degrees(angle)
