import numpy as np

class LP:
    def __init__(self, A, b, c):
        self.A = np.array(A, dtype=float)
        self.b = np.array(b, dtype=float)
        self.c = np.array(c, dtype=float)
        self.m, self.n = self.A.shape

    def is_valid_basis(self, basis):
        if len(basis) != self.m:
            return False
        B = self.A[:, basis]
        return np.linalg.matrix_rank(B) == self.m


# read_init() forms the original LP based on user input
def read_init() :
    print('An LP is of the form max{cTx: Ax<=b, x>=0}')
    n = int(input('Number of variables:'))
    m = int(input('Number of constraints:'))
    A = np.zeros((m, n))
    b = np.zeros(m)
    for i in range(m):
        while True:
            try:
                row = list(map(float, input(f"Enter coefficients for row {i+1} of A:").split()))
                if len(row) != n:
                    print(f"Expected {n} values, got {len(row)}. Try again.")
                    continue
                break
            except ValueError:
                print(f"Invalid input, enter {n} values seperated by spaces")
                print("Example: 1 2 3 4")
        A[i] = row
        cnst = float(input(f"Enter value for b_{i+1}:"))
        b[i] = cnst
    
    while True:
        try:
            c = list(map(float, input(f"Enter all {n} costs for c:").split()))
            if len(c) != n:
                print(f"Expected {n} values, got {len(row)}. Try again.")
                continue
            break
        except ValueError:
            print(f"Invalid input, enter {n} values seperated by spaces")
            print("Example: 1 2 3 4")
    
    lp = LP(A, b, c)
    return lp

# fmt(x) returns x as a coefficient for a linear equation. Returns an empty string if x = 1.
def fmt(x):
    if x == 1:
        return ""
    elif x == -1:
        return "-"
    else:
        return f"{x:g}"

# print_LP(lp) prints an LP as a system of linear inequalities.
def print_LP(lp):
    print("LP: max " + " + ".join(fmt(lp.c[i]) + f"x{i+1}" for i in range(lp.n) if not np.isclose(lp.c[i], 0)))
    print("s.t.")
    for j in range(lp.m):
        terms = []
        for i in range(lp.n):
            if lp.A[j, i] == 0:
                continue
            if terms and lp.A[j, i] < 0:
                terms.append("- " + fmt(abs(lp.A[j, i])) + f"x{i+1}")
            else:
                terms.append(fmt(lp.A[j, i]) + f"x{i+1}")
        print(" + ".join(terms).replace("+ -", "-") + f" <= {lp.b[j]:g}")