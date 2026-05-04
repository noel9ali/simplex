from lp import LP, read_init
import numpy as np

# GLOBAL CONSTANTS
DEGENERATE = 300
UNBOUNDED  = 400
SOLVED     = 500
INFEASIBLE = 600


# iterate_basis(basis, lp) finds an entering and leaving variable for an LP
# in canonical form using Bland's rule
def iterate_basis(basis, lp, verbose=True):
    # find an entering variable
    k = -1
    for i in range(lp.n):
        if lp.c[i] > 0:
            k = i
            break
    # if no entering variable found, return original basis
    if k == -1:
        if verbose:
            print("No entering variables exist, basis is optimal.")
        return SOLVED

    ratios = np.full(lp.m, np.inf)
    for j in range(lp.m):
        if lp.A[j, k] <= 0:
            continue
        else:
            ratios[j] = lp.b[j] / lp.A[j, k]

    l = np.argmin(ratios)
    if np.isinf(ratios[l]):
        return UNBOUNDED
    else:
        if verbose:
            print(f"Leaving variable: {basis[l]}")
            print(f"Entering variable: {k}")
            print(f"New basis: " + ", ".join(f"{basis[i]}" for i in range(lp.m)))
        basis = [x for x in basis if x != basis[l]] + [k]
        basis.sort()
        return basis


# simplex(basis, lp) performs simplex iterations on lp until a result is found.
#   If optimal, returns the optimal basis; otherwise returns UNBOUNDED.
#   verbose=False suppresses all output (used during Phase I).
def simplex(basis, lp, verbose=True):
    if verbose:
        print("Performing simplex with starting basis: " + ", ".join(f"{basis[i]}" for i in range(lp.m)))
    while True:
        canonical = lp.canonical_form(basis)
        if verbose:
            print("The canonical form is:")
            canonical.display()
        result = iterate_basis(basis, canonical, verbose=verbose)

        if result == SOLVED:
            return basis
        elif result == UNBOUNDED:
            return UNBOUNDED
        else:
            basis = result


# print_result(basis, lp) runs simplex and prints the optimal solution, or
#   reports that the LP is unbounded.
def print_result(basis, lp):
    result = simplex(basis, lp)
    if result == UNBOUNDED:
        print("The LP is unbounded.")
    else:
        x_bar = "x\u0305"
        x, obj = lp.get_solution(result)
        print(f"Optimal solution: {x_bar} = (" + ", ".join(f"{x[i]:g}" for i in range(lp.n)) + ")")
        print(f"Optimal value: {obj:g}")


# initial_feasible(lp) finds an initial feasible basis for lp via the
#   auxiliary problem method (Big-M / two-phase setup)
def initial_feasible(lp):

    m, n = lp.A.shape
    aux = lp.SEF().aux()
    basis = list(range(n + m, n + 2*m))
    result = simplex(basis, aux, verbose=False)
    if result == INFEASIBLE or result == UNBOUNDED:
        return INFEASIBLE
    else:
        x, obj = aux.get_solution(result)
    if not np.isclose(obj, 0):
        return INFEASIBLE
    
    clean_basis = [b for b in result if b < n + m]

    # if basis is short, find replacement columns
    for j in range(n + m):
        if len(clean_basis) == m:
            break
        if j not in clean_basis:
            clean_basis.append(j)

    clean_basis.sort()
    return clean_basis


# solve(lp) runs the full two-phase simplex on lp and prints the result
#   in terms of the original variables (excluding SEF slack variables).
def solve(lp):
    n = lp.n  # original variable count, before SEF adds slacks
    basis = initial_feasible(lp)
    if basis == INFEASIBLE:
        print("The LP is infeasible.")
        return
    sef = lp.SEF()
    result = simplex(basis, sef)
    if result == UNBOUNDED:
        print("The LP is unbounded.")
        return
    x, obj = sef.get_solution(result)
    x_bar = "x\u0305"
    print(f"Optimal solution: {x_bar} = (" + ", ".join(f"{x[i]:g}" for i in range(n)) + ")")
    print(f"Optimal value: {obj:g}")
