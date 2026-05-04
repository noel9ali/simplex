import numpy as np


# fmt(x) returns x as a coefficient for a linear equation. Returns an empty string if x = 1.
def fmt(x):
    if np.isclose(x, 1):
        return ""
    elif np.isclose(x, -1):
        return "-"
    else:
        return f"{x:g}"


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

    # SEF() translates the LP into its standard equality form
    def SEF(self):
        m, n = self.A.shape
        A_sef = np.zeros((m, n + m))
        A_sef[:, :n] = self.A
        A_sef[:, n:] = np.eye(m)
        c_sef = np.concatenate([self.c, np.zeros(m)])
        b_sef = self.b.copy()
        for i in range(m):
            if self.b[i] < 0:
                b_sef[i] *= -1
                A_sef[i] *= -1
        return LP(A_sef, b_sef, c_sef)

    # aux() adds auxiliary variables to the LP (assumed to be in SEF)
    def aux(self):
        m, n = self.A.shape
        A_aux = np.hstack([self.A, np.eye(m)])
        c_aux = np.concatenate([np.zeros(n), -np.ones(m)])
        aux = LP(A_aux, self.b, c_aux)
        return aux

    # display() prints the LP as a system of linear inequalities
    def display(self):
        print("LP: max " + " + ".join(fmt(self.c[i]) + f"x{i+1}" for i in range(self.n) if not np.isclose(self.c[i], 0)))
        print("s.t.")
        for j in range(self.m):
            terms = []
            for i in range(self.n):
                if np.isclose(self.A[j, i], 0):
                    continue
                if terms and self.A[j, i] < 0:
                    terms.append("- " + fmt(abs(self.A[j, i])) + f"x{i+1}")
                else:
                    terms.append(fmt(self.A[j, i]) + f"x{i+1}")
            print(" + ".join(terms).replace("+ -", "-") + f" <= {self.b[j]:g}")

    # canonical_form(basis) returns the canonical form of the LP for the given basis
    def canonical_form(self, basis):
        # check if valid basis
        if not self.is_valid_basis(basis):
            print(" ".join(f"{basis[i]}" for i in range(len(basis))) + " is not a valid basis.")
            return None

        # find inverse of B
        B = self.A[:, basis]
        B_inv = np.linalg.inv(B)

        # find canonical A, b, c
        A_prime = B_inv @ self.A
        b_prime = B_inv @ self.b
        c_basic = self.c[basis]
        c_prime = self.c - c_basic @ A_prime

        canonical_lp = LP(A_prime, b_prime, c_prime)
        return canonical_lp

    # get_solution(basis) returns the primal solution and objective value for the given basis
    def get_solution(self, basis):
        B = self.A[:, basis]
        x = np.zeros(self.n)
        x[basis] = np.linalg.solve(B, self.b)
        obj = self.c @ x
        return x, obj


# read_init() forms an LP from user input
def read_init():
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
