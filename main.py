from lp import read_init
from simplex import solve

if __name__ == "__main__":
    lp = read_init()
    solve(lp)