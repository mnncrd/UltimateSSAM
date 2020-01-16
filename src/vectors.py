"""Operations on vectors module.

This module performs operations on vectors. For instance, it computes the vector
norm, the scalar and cross products, as well as sum and difference between
vectors.
"""

# Import libraries
import math

def vector_norm(v_i):

    """Computes the vector norm.

    Computes the norm of the given vector.

    Args:
        v_i: A vector of size 3.

    Returns:
        A float equalling the vector norm.
    """

    v_i_norm = math.sqrt(v_i[0]**2 + v_i[1]**2 + v_i[2]**2)
    return v_i_norm

def scal_product(v_i, v_j):

    """Computes the scalar product.

    Computes the scalar product, or dot product, between the two given vectors.

    Args:
        v_i: A vector of size 3.
        v_j: A vector of size 3.

    Returns:
        A float equalling the scalar product.
    """

    scal_prod = v_i[0]*v_j[0] + v_i[1]*v_j[1] + v_i[2]*v_j[2]
    return scal_prod

def cross_product(v_i, v_j):

    """Computes the cross product.

    Computes the cross product, or vector product, between the two given
    vectors.

    Args:
        v_i: A vector of size 3.
        v_j: A vector of size 3.

    Returns:
        A vector equalling the cross product.
    """

    v_k_x = v_i[1]*v_j[2]-v_i[2]*v_j[1]
    v_k_y = v_i[2]*v_j[0]-v_i[0]*v_j[2]
    v_k_z = v_i[0]*v_j[1]-v_i[1]*v_j[0]
    cross_prod = [v_k_x, v_k_y, v_k_z]
    return cross_prod

def compute_frac_vect(v_i, div):

    """Computes the quotient of a vector by a divisor.

    Divides the given vector by the given divisor.

    Args:
        v_i: A vector of size 3.
        div: A float by which to divide the vector.

    Returns:
        A vector equalling the quotient.
    """

    v_div = [i/div for i in v_i]
    return v_div

def compute_diff_vect(v_i, v_j):

    """Computes the diffence between two vectors.

    Substracts each dimension of the second given vector to the dimension of
    the first one.

    Args:
        v_i: A vector of size 3.
        v_j: A vector of size 3.

    Returns:
        A vector equalling the difference.
    """

    v_ij = [i-j for i, j in zip(v_i, v_j)]
    return v_ij

def compute_sum_vect(v_i, v_j):

    """Computes the sum between two vectors.

    Adds each dimension of the second given vector to the dimension of the first
    one.

    Args:
        v_i: A vector of size 3.
        v_j: A vector of size 3.

    Returns:
        A vector equalling the sum.
    """

    v_ij = [i+j for i, j in zip(v_i, v_j)]
    return v_ij
