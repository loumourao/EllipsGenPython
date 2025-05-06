# TODO: Either extinguish or adapt the varargout stuff to Python as I'm not sure
#       whether scipy root_scalar return values match those of MATLAB's 'fzeros'.
#       Adapt documentation to Python's standards
#       TEST

import numpy as np

from matrices import matrices
from numpy import linalg as LA
from scipy.optimize import root_scalar

# Computes the contact points between two ellipses.
#
# Given two ellipses defined by their parameter sets E_i and E_j, this 
# function finds the contact points x_i and x_j that satisfy the geometric 
# conditions.
#
# -------------------------------------------------------------------------
# Inputs:
#   E_i - Parameter set of the first ellipse [a, b, theta, o_x, o_y]
#   E_j - Parameter set of the second ellipse [a, b, theta, o_x, o_y]
#
# Outputs:
#   x_i - Contact point on the first ellipse
#   x_j - Contact point on the second ellipse
#   varargout:
#     {1} - Residual error for the first contact equation (fval1)
#     {2} - Residual error for the second contact equation (fval2)
#     {3} - Exit flag for the first solver (exitflag1)
#     {4} - Exit flag for the second solver (exitflag2)
#     {5} - Solver output details for the first contact point (output1)
#     {6} - Solver output details for the second contact point (output2)
#
# -------------------------------------------------------------------------

def contact_solve(E_i, E_j):
    # Compute quadratic matrices and centers of the ellipses
    Q_i, o_i, *_ = matrices(E_i)
    Q_j, o_j, *_ = matrices(E_j)

    # Define quadratic potential functions for each ellipse
    ei = lambda x : x.conj().T @ Q_i @ x - 1;  # Ellipse i equation
    ej = lambda x : x.conj().T @ Q_j @ x - 1;  # Ellipse j equation

    # Define the parametric form w(u) of co-gradient locus H_ij
    # WARNING: THIS ONE MUST BE TESTED TO ENSURE CORRECT VALUES
    #          IT MAY BE LEAST SQUARES INSTEAD
    #          solve - for square and non-singular A matrix
    #          lstsq - when A is not square 
    w = lambda u : LA.solve((1 - u) @ Q_i + u @ Q_j,
                            (1 - u) @ Q_i @ o_i + u @ Q_j @ o_j)

    # Define contact conditions for each ellipse
    contact_point_1 = lambda u: ei(w(u) - o_i)
    contact_point_2 = lambda u: ej(w(u) - o_j)

    parameter_1 = root_scalar(contact_point_1, bracket=[0, 1], method='brenqt').root
    x_i = w(parameter_1)

    parameter_2 = root_scalar(contact_point_2, bracket=[0, 1], method='brenqt').root
    x_j = w(parameter_2)

    return x_i, x_j
