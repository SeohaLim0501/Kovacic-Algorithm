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
