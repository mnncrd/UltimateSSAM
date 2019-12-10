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