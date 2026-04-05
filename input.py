import numpy as np

class LP:
    def __init__(self, A, b, c):
        self.A = np.array(A, dtype=float)
        self.b = np.array(b, dtype=float)
        self.c = np.array(c, dtype=float)
        self.m, self.n = self.A.shape

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
    return "" if x == 1 else f"{x:g}"

# print_LP(lp) prints an LP as a system of linear inequalities.
def print_LP(lp):
    print("LP: max " + " + ".join(fmt(lp.c[i]) + f"x{i+1}" for i in range(lp.n) if not lp.c[i] == 0))
    print("s.t.")
    for j in range(lp.m):
        print(" + ".join(fmt(lp.A[j, i]) + f"x{i+1}" for i in range(lp.n) if not lp.A[j, i] == 0) + f" <= {lp.b[j]:g}")