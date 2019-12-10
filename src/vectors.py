"""Operations on vector module.

This module performs operations on vectors. For instance, it computes the
vector norm.
"""

import math

def vector_norm(v_i):

    """Computes the vector norm"""

    v_i_norm = math.sqrt(v_i[0]**2 + v_i[1]**2 + v_i[2]**2)
    return v_i_norm

def scal_product(v_i, v_j):

    """Computes the scalar product between two vectors"""

    scal_prod = v_i[0]*v_j[0] + v_i[1]*v_j[1] + v_i[2]*v_j[2]
    return scal_prod

def compute_diff_vect(v_i, v_j):

    """Computes the diffence between two vectors"""

    v_ij = [i-j for i, j in zip(v_i, v_j)]
    return v_ij
