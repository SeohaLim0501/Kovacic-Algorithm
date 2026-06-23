# Kovacic-Algorithm
My personal project to perform Kovacic Algorithm.

## Setup

```bash
git clone https://github.com/SeohaLim0501/Kovacic-Algorithm.git
cd Kovacic-Algorithm
```

Create and activate a Miniconda environment:

```bash
conda create -n kovacic python=3.11
conda activate kovacic
```

Install dependencies:

```bash
pip install -r requirements.txt
```

## Run

### Default Examples

```bash
python main.py
```

By default, `main.py` reads test cases from `examples.txt`.

### Integration Timeout

Integration timeout defaults to 1 second. Use `-t` or `--timeout` to change it:

```bash
python main.py -t 5
```

### Equation Solve Timeout

The timeout for solving the Case 3 omega equation defaults to 1 second. Use
`-t2` or `--timer2` to change it:

```bash
python main.py -t2 20
```

Both timeout options can be used together:

```bash
python main.py -t 5 -t2 20
```

### Input Text Files

Use `-c` or `--cases` to choose a different text file:

```bash
python main.py -c input.txt
```

For your own input, you can either write to `input.txt` from the shell or create/edit any `.txt` file yourself. Then pass that file with `-c`:

```bash
python main.py -c input.txt
```

Each non-comment line in the input file can be either:

```txt
x + 77/(16*x**2)
My test | (4*x**2 + 8*x + 6)/(2*x + 1)**2
```

Second-order ODE input is also accepted. Use SymPy derivative syntax with
`y(x)` or `z(x)`; both `Eq` and left-hand-side-only forms are supported:

```txt
ODE with Eq | Eq(diff(y(x), x, 2) + 2*x*diff(y(x), x) + y(x), 0)
ODE lhs only | diff(y(x), x, 2) - y(x)
```

ODE input is regularized automatically and solved through `SecOrdSolver`.

In PowerShell, overwrite `input.txt` with one test case:

```powershell
"x + 77/(16*x**2)" | Set-Content input.txt
```

Or include a title:

```powershell
"My test | x + 77/(16*x**2)" | Set-Content input.txt
```

Append another test case:

```powershell
"(4*x**2 + 8*x + 6)/(2*x + 1)**2" | Add-Content input.txt
```

Clear all contents of `input.txt`:

```powershell
Clear-Content input.txt
```

## Testing with `test.py`

`test.py` runs the equations in `kovacic_test.txt` and checks whether the
detected Kovacic case matches the expected case stored on each line. It also:

- compiles the project modules before running the tests;
- continues to the next test when one test raises an exception;
- prints green `Correct!` or red `Incorrect!` for each result;
- reports the score and the numbers of all incorrect tests;
- reports the average elapsed time for each expected case;
- writes the complete output automatically to `result.txt` without ANSI
  color codes.

Run the default test suite with:

```bash
python test.py
```

The integration and Case 3 equation-solve timeouts use the same options as
`main.py`:

```bash
python test.py -t 5 -t2 20
```

Use another test file or output file with `-c` and `-o`:

```bash
python test.py -c my_tests.txt -o my_results.txt
```

`result.txt` is overwritten on each default run.

### Editing `kovacic_test.txt`

Each non-comment line must contain exactly three fields separated by `|`:

```text
title | r(x) | expected_case
```

For example:

```text
Constant Case 1 | 1 | 1
Airy-like Case 4 | -x**3 + 3*x - 1 | 4
Trivial equation | 0 | 0
```

The expected case must be an integer from `0` to `4`:

- `0`: the trivial equation `z'' = 0`;
- `1`: Kovacic Case 1;
- `2`: Kovacic Case 2;
- `3`: Kovacic Case 3;
- `4`: Cases 1, 2, and 3 all fail, so no Liouvillian solution is found.

Blank lines and lines beginning with `#` are ignored:

```text
# title | r(x) | expected_case
[case 2] example | x + 77/(16*x**2) | 2
```

The expression field is parsed as a SymPy expression. Use `I` for the
imaginary unit, `S(1)/2` for an explicit SymPy rational when needed, and `x`
as the independent variable. General second-order ODE expressions may also
use `Eq`, `diff`, `Derivative`, `y(x)`, and `z(x)`.

## `kovacic.py`

The `kovacic.py` module exposes `simpleKovacic`, which applies Kovacic's
algorithm directly to an equation already written in reduced form:

```text
z'' = r(x)z
```

Import it with:

```python
from sympy import symbols
from kovacic import simpleKovacic

x = symbols("x")
r = (x**2 + 1) / x**2
result = simpleKovacic(r, x)
```

### `simpleKovacic` Input

The function signature is:

```python
simpleKovacic(
    r,
    x,
    integration_timeout_seconds=1,
    solve_timeout_seconds=1,
)
```

- `r` is the SymPy rational function on the right-hand side of
  `z'' = r(x)z`.
- `x` is the SymPy `Symbol` used as the independent variable.
- `integration_timeout_seconds` limits the integration of `omega` in Cases 1
  and 2. If integration times out, the solution remains in formal
  `exp(Integral(...))` form.
- `solve_timeout_seconds` limits solving the algebraic equation for `omega`
  in Case 3.

`simpleKovacic` expects an already reduced equation. To start with a general
homogeneous second-order equation such as `y'' + a*y' + b*y = 0`, use
`SecOrdSolver` from `SecOrdSolver.py`.

### `simpleKovacic` Output

The return value is a dictionary with three top-level keys:

```python
result = {
    "solution": ...,  # solution expression, Case 3 data, or failure message
    "status": ...,    # short description of the outcome
    "debug": {...},   # pole analysis and case-specific intermediate data
}
```

Not every key inside `debug` is always present. Case-specific data is added
only when that case is valid and reached, and the function returns as soon as
a case succeeds.

#### `solution`

For the trivial equation `z'' = 0`, `solution` is the general linear
solution:

```python
C1*x + C2
```

For a successful Case 1 or Case 2, it is a SymPy expression representing one
Liouvillian solution. When the final integration times out, the expression
may contain an unevaluated `Integral`.

For a successful Case 3, a closed-form root for `omega` may not be available,
so `solution` has this structure:

```python
{
    "form": "exp(Integral(omega, x))",
    "omega_symbol": omega,
    "omega_equation": Eq(..., 0),
}
```

If Cases 1, 2, and 3 all fail, the value is:

```python
"No Liuvillian solution for this!"
```

#### `status`

`status` is a human-readable summary. Current values include:

```text
Solved trivial equation z'' = 0
Solved by Kovacic case 1
Solved by Kovacic case 2
Solved by Kovacic case 3 with n=4
Solved by Kovacic case 3 with n=6
Solved by Kovacic case 3 with n=12
No Liuvillian solution found
```

#### `debug`

For every nontrivial input, `debug` begins with a shared pole analysis:

```python
debug["pole_analysis"] = {
    "num": ...,
    "den": ...,
    "deg_num": ...,
    "deg_den": ...,
    "finite_poles": {pole: order, ...},
    "num_finite_poles": ...,
    "pole_at_infinity": ...,
    "case_analysis": {
        "case1_valid": True_or_False,
        "case2_valid": True_or_False,
        "case3_valid": True_or_False,
    },
}
```

The remaining keys depend on the Kovacic case.

##### Case 1 debug data

When Case 1 is valid, the finite-pole and infinity data are stored as:

```python
debug["case1_finite_data"] = {
    "sqrtR": {pole: value, ...},
    "alphaPlus": {pole: value, ...},
    "alphaMinus": {pole: value, ...},
}

debug["case1_infinity_data"] = {
    "sqrtR_inf": ...,
    "alphaPlus_inf": ...,
    "alphaMinus_inf": ...,
}
```

`debug["case1_candidates"]` is a list of candidate dictionaries:

```python
{
    "d": ...,
    "d_expr": ...,
    "omega": ...,
    "finite_signs": {pole: 1_or_minus_1, ...},
    "inf_sign": 1_or_minus_1,
    "finite_alpha": {pole: value, ...},
    "infinite_alpha": ...,
}
```

On success, the following keys are also present:

```python
debug["successful_case"] = 1
debug["successful_candidate"] = {...}
debug["case1_step3"] = {
    "p": ...,
    "omega": ...,
    "z": ...,
    "z_evaluated": ...,
    "integration_timed_out": True_or_False,
    "coeff_solution": {...},
}
```

If the case does not produce a candidate, `debug["case1_status"]` explains
why. In the current implementation, the `d == 0` success path contains the
misspelled key `coeff_soution` instead of `coeff_solution`.

##### Case 2 debug data

When Case 2 is valid, its pole data is stored as:

```python
debug["case2_finite_data"] = {
    "E": {pole: [values], ...},
    "bval": {pole: value, ...},
}

debug["case2_infinity_data"] = {
    "E_inf": [values],
    "b_inf": ...,
}
```

`debug["case2_candidates"]` contains dictionaries of the form:

```python
{
    "d": ...,
    "d_expr": ...,
    "theta": ...,
    "e_finite": {pole: value, ...},
    "e_inf": ...,
}
```

On success, the additional data is:

```python
debug["successful_case"] = 2
debug["successful_candidate"] = {...}
debug["case2_step3"] = {
    "p": ...,
    "theta": ...,
    "phi": ...,
    "omega": ...,
    "z": ...,
    "z_evaluated": ...,
    "integration_timed_out": True_or_False,
    "coeff_solution": {...},
    "candidate": {...},
}
```

If Case 2 is invalid or produces no candidate, `debug["case2_status"]`
contains the reason.

##### Case 3 debug data

Case 3 tests `n = 4`, `6`, and `12`. Its finite-pole and infinity data are
therefore dictionaries keyed by `n`:

```python
debug["case3_finite_data"] = {
    n: {
        "E": {pole: [values], ...},
        "bval": {pole: value, ...},
    },
    ...,
}

debug["case3_infinity_data"] = {
    n: {
        "E_inf": [values],
        "b_inf": ...,
    },
    ...,
}
```

The values of `n` are evaluated lazily. If Case 3 succeeds for `n = 4`, the
dictionaries contain only the key `4`; keys `6` and `12` are added only when
those values are actually reached.

`debug["case3_candidates"]` contains dictionaries of the form:

```python
{
    "n": 4_or_6_or_12,
    "d": ...,
    "d_expr": ...,
    "theta": ...,
    "S": ...,
    "e_finite": {pole: value, ...},
    "e_inf": ...,
}
```

When a recurrence polynomial is found, Case 3 adds:

```python
debug["case3_N"] = ...
debug["case3_Efinite"] = {...}
debug["case3_Einfinite"] = ...
debug["case3_P0"] = ...
debug["case3_omega_equation"] = ...
debug["case3_solve_timed_out"] = True_or_False
debug["case3_omega_solutions"] = [...]  # present when solving succeeds

debug["successful_case"] = 3
debug["successful_candidate"] = {...}
debug["case3_step3"] = {
    "P": ...,
    "theta": ...,
    "S": ...,
    "n": ...,
    "d": ...,
    "P_sequence": {degree: expression, ...},
    "coeff_solution": {...},
    "candidate": {...},
}
```

`debug["case3_status"]` reports whether the case was invalid, had no
nonnegative integer `d`, failed to find a recurrence polynomial, timed out
while solving for `omega`, or solved the `omega` equation.

## Output Format

The program first prints the selected input file and timeout settings:

```text
========================================================================
Simple Kovacic Algorithm
Input file: examples.txt
Integration timeout: 1 second
Equation Solve timeout: 1 second
========================================================================
```

Each test result then contains:

```text
------------------------------------------------------------------------
Test 1: Example title
Elapsed: 0.123 seconds
r(x): ...
Status: Solved by Kovacic case 1
Solution: ...
Pole analysis:
  finite poles: 0 (order 2)
  pole at infinity: order 1
Valid cases: case 1: yes, case 2: no, case 3: no
```

- `Elapsed` is the total running time for that test case.
- `r(x)` is the normalized right-hand side of `z'' = r(x)z`.
- `Status` reports whether a Kovacic case found a solution.
- `Solution` is the resulting expression, or `not found`.
- `finite poles` lists each finite pole and its order.
- `pole at infinity` gives the order used by the case analysis.
- `Valid cases` shows which Kovacic cases satisfy the pole conditions.

When a case succeeds, an additional case-specific section is printed.

### Case 1 Details

Case 1 prints the selected alpha value at each finite pole and at infinity,
followed by `d` and `omega`:

```text
Case 1 details:
  alpha values:
    finite_alpha: 0: -3/2
    infinite_alpha: 1/2
  d: 2
  omega: ...
```

### Case 2 Details

Case 2 prints the finite and infinite E sets, followed by `d` and `omega`:

```text
Case 2 details:
  E sets:
    E_0 = [1, 2, 3]
    E_infinity = [1]
  d: 0
  omega: ...
```

If there are several finite poles, one line is printed for each set, such as
`E_0`, `E_1`, and so on.

### Case 3 Details

Case 3 prints the E sets selected for its value of `N`, followed by `d`, `N`,
and the recurrence polynomial `P_0`:

A Case 3 solution usually is not returned directly in closed form. Instead,
the program provides an algebraic equation for `omega`. For any root `omega`
of that equation, a corresponding solution can be constructed as
`z = exp(Integral(omega, x))`.

```text
Case 3 details:
  E sets:
    E_0 = [1, 2, 3]
    E_infinity = [4, 5]
  d: 2
  N: 12
  P_0: ...
```
