"""Operations on vector module.

This module performs operations on vectors. For instance, it computes the
vector norm.
"""

import math

def vector_norm(v_i):

    """Computes the vector norm"""

    v_i_norm = math.sqrt(v_i[0]**2 + v_i[1]**2 + v_i[2]**2)
    return v_i_norm