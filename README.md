# Simplex LP Solver

A two-phase Simplex Method implementation for solving linear programs in Python.

## Problem Form

The solver accepts LPs of the form:

```
max { cᵀx : Ax ≤ b, x ≥ 0 }
```

Internally, the LP is converted to Standard Equality Form (SEF) by adding slack variables before solving.

## Usage

Run from the command line:

```bash
python main.py
```

You will be prompted to enter:
- Number of variables `n` and constraints `m`
- Each row of the constraint matrix `A`
- Each right-hand side value `b`
- The cost vector `c`

### Example

For the LP `max 3x₁ + 2x₂` subject to `x₁ + x₂ ≤ 4`, `2x₁ + x₂ ≤ 6`:

```
Number of variables: 2
Number of constraints: 2
Enter coefficients for row 1 of A: 1 1
Enter value for b_1: 4
Enter coefficients for row 2 of A: 2 1
Enter value for b_2: 6
Enter all 2 costs for c: 3 2
```

Output:
```
Optimal solution: x̅ = (2, 2)
Optimal value: 10
```

## Project Structure

| File | Purpose |
|---|---|
| `lp.py` | `LP` class (data model, SEF/aux transforms, canonical form) and `read_init()` |
| `simplex.py` | Simplex algorithm: Phase I (`initial_feasible`), Phase II (`simplex`, `iterate_basis`), and top-level `solve` |
| `main.py` | Entry point |

## Requirements

- Python 3
- NumPy

Install NumPy with:

```bash
pip install numpy
```
