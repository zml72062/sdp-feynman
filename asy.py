import sympy
import scipy.spatial as sp
from typing import Tuple, Any


def dimred(points: sympy.Matrix) -> Tuple[sympy.Matrix, sympy.Matrix, sympy.Matrix]:
    r"""
    Given a set of points [v1, ..., vn] in k-dimensional space, find a 
    superplane that contains all n points yet have the minimum dimension.
    
    Arguments:
        points (sympy.Matrix) - a (k*n) matrix whose i-th column is the
        coordinate of the i-th point
    
    Returns:
        reduced_points (sympy.Matrix) - a (r*n) matrix whose i-th column
        is the reduced coordinate of the i-th point. The dimension r equals
        the dimension of the superplane

        transform_weight (sympy.Matrix) - a (k*r) column orthogonal matrix
        describing the coordinate frame on the r-dimensional superplane

        transform_bias (sympy.Matrix) - a (k*1) matrix describing how to
        move the k-dimensional origin to the origin on the superplane

    Note: To recover the i-th input point, one only need to do

        transform_weight * reduced_points[:, i] + transform_bias
    """
    transform_bias = points[:, 0]
    for i in range(points.cols):
        points[:, i] -= transform_bias
    transform_weight, reduced_points = points.QRdecomposition()
    return reduced_points, transform_weight, transform_bias


def getvecs(points: sympy.Matrix) -> sympy.Matrix:
    r"""
    Get the scaling vectors (r_1, r_2, ..., r_n) corresponding to the point
    set determined by terms in polynomial UF. 

    To be specific, each row of the result matrix is a scaling vector (r_1, 
    r_2, ..., r_n) of Feynman parameters, which corresponds to (x_1 ~ 
    \rho^{r_1}, ..., x_n ~ \rho^{r_n}), where \rho is the infinitesimal 
    scaling parameter. The returned scaling vectors are those that integrate
    to non-zero leading contributions.
    """
    reduced_points, _, _ = dimred(points)
    hull = sp.ConvexHull(reduced_points.T)
    vectors = set()
    syms = sympy.symbols(",".join([f"x{i}" for i in range(points.rows)]))
    for simplex in hull.simplices:
        selected_points = points[:, simplex.tolist()].T.copy()
        bias = -selected_points[:, -1]
        for i in range(selected_points.rows):
            selected_points[i, -1] = -1
        solutions = sympy.linsolve((selected_points, bias), syms)
        for s in solutions: # check direction
            s = s.subs([(sym, 0) for sym in syms])
            orders = points[:-1, :].T.multiply(sympy.Matrix(s[:-1])).T + points[-1, :]
            if all(t >= s[-1] for t in orders):
                vectors.add(s[:-1])
    
    return sympy.Matrix(list(vectors))
        
        
def run(point_list: Any):
    return getvecs(sympy.Matrix(point_list).T)

