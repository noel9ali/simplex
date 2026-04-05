from input import LP, print_LP
import numpy as np

# GLOBAL CONSTANTS
DEGENERATE=300
UNBOUNDED=400
SOLVED=500

# find_canonical(B, lp) returns the canonical form of an LP for the basis B.
def find_canonical(basis, lp):
    # check if valid basis
    if not lp.is_valid_basis(basis):
        print(" ".join(f"{basis[i]}" for i in range(len(basis))) + " is not a valid basis.")
        return None

    # find inverse of B
    B = lp.A[:, basis]
    B_inv = np.linalg.inv(B)

    # find canonical A, b, c
    A_prime = B_inv @ lp.A
    b_prime = B_inv @ lp.b
    c_basic = lp.c[basis]
    c_prime = lp.c - c_basic @ A_prime

    canonical_lp = LP(A_prime, b_prime, c_prime)
    print("The canonical form is:")
    print_LP(canonical_lp)
    return canonical_lp

# iterate_basis(basis, lp) finds an entering and leaving variable for an LP in canonical form using Bland's rule
def iterate_basis(basis, lp):
    # find an entering variable
    k = -1
    for i in range(lp.n):
        if lp.c[i] > 0:
            k = i
            break
    # if no entering variable found, return original basis
    if k == -1:
        print("No entering variables exist, basis is optimal.")
        return SOLVED
    
    ratios = np.full(lp.m, np.inf)
    for j in range(lp.m):
        if lp.A[j, k] <= 0:
            continue
        else:
            ratios[j] = lp.b[j] / lp.A[j, k]
    
    l = np.argmin(ratios)
    if ratios[l] == np.inf:
        return UNBOUNDED
    else:
        print(f"Leaving variable: {l}")
        print(f"Entering variable: {k}")
        print(f"New basis: " + ", ".join(f"{basis[i]}" for i in range(lp.m)))
        basis = [x for x in basis if x != basis[l]] + [k]
        basis.sort()
        return basis

# simplex(basis, lp) performs simplex iterations on lp until a result is found. 
#   If optimal, it returns the optimal basis
def simplex(basis, lp):
    print("Performing simplex with starting basis: " + ", ".join(f"{basis[i]}" for i in range(lp.m)))
    while True:
        canonical = find_canonical(basis, lp)
        result = iterate_basis(basis, canonical)
        
        if result == SOLVED:
            return basis
        elif result == UNBOUNDED:
            return UNBOUNDED
        else:
            basis = result

def get_solution(basis, lp):
    B = lp.A[:, basis]
    x = np.zeros(lp.n)
    x[basis] = np.linalg.solve(B, lp.b)
    obj = lp.c @ x
    return x, obj

def print_result(basis, lp):
    result = simplex(basis, lp)
    if result == UNBOUNDED:
        print("The LP is unbounded.")
    else:
        x_bar = "x\u0305"
        x, obj = get_solution(result, lp)
        print(f"Optimal solution: {x_bar} = (" + ", ".join(f"{x[i]:g}" for i in range(lp.n)) + ")")
        print(f"Optimal value: {obj:g}")