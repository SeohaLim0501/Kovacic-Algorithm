"""Utilities for handling second-order differential equations."""

from sympy import (
    Derivative,
    Eq,
    Function,
    cancel,
    diff,
    expand,
    simplify,
    symbols,
    together,
)
from sympy.core.expr import Expr
from sympy.core.function import AppliedUndef
from sympy.core.sympify import SympifyError, sympify
from kovacic import simpleKovacic


def regularize(expr):
    """Return ``r`` after reducing ``y'' + a*y' + b*y = 0`` to ``z'' = r*z``.

    The old and new dependent variables are related by
    ``y = exp(-Integral(a/2, x))*z``, and the returned rational function is
    ``r = diff(a, x)/2 + a**2/4 - b``.

    ``expr`` may be an ``Eq`` or a left-hand-side expression understood to
    equal zero.  A ``ValueError`` is raised when the input is not a
    homogeneous linear second-order equation with rational coefficients.
    """
    try:
        equation = sympify(expr)
    except (SympifyError, TypeError, ValueError) as error:
        raise ValueError("Could not interpret the differential equation.") from error

    if isinstance(equation, Eq):
        equation = equation.lhs - equation.rhs
    elif not isinstance(equation, Expr):
        raise ValueError("Input must be a SymPy Eq or expression.")

    equation = expand(equation)
    second_derivatives = {
        derivative
        for derivative in equation.atoms(Derivative)
        if len(derivative.variables) == 2
        and derivative.variables[0] == derivative.variables[1]
        and isinstance(derivative.expr, AppliedUndef)
        and len(derivative.expr.args) == 1
        and derivative.expr.args[0] == derivative.variables[0]
    }
    if len(second_derivatives) != 1:
        raise ValueError("The equation must contain exactly one second derivative.")

    y2 = second_derivatives.pop()
    y = y2.expr
    x = y2.variables[0]
    y1 = diff(y, x)

    if equation.atoms(AppliedUndef) != {y}:
        raise ValueError("The equation must contain exactly one unknown function.")
    if not equation.atoms(Derivative) <= {y1, y2}:
        raise ValueError("Only the first and second derivatives are allowed.")

    coefficient_y2 = equation.coeff(y2)
    coefficient_y1 = equation.coeff(y1)
    without_derivatives = simplify(
        equation - coefficient_y2*y2 - coefficient_y1*y1
    )
    coefficient_y = without_derivatives.coeff(y)
    remainder = simplify(without_derivatives - coefficient_y*y)

    coefficients = (coefficient_y2, coefficient_y1, coefficient_y)
    if coefficient_y2 == 0 or any(c.has(y, y1, y2) for c in coefficients):
        raise ValueError("The equation must be linear in the unknown function.")
    if remainder != 0:
        raise ValueError("The equation must be homogeneous.")
    if not all(c.is_rational_function(x) is True for c in coefficients):
        raise ValueError("All coefficients must be rational functions.")

    a = cancel(together(coefficient_y1 / coefficient_y2))
    b = cancel(together(coefficient_y / coefficient_y2))
    r = cancel(together(diff(a, x)/2 + a**2/4 - b))
    return r


def SecOrdSolver(
    expr,
    integration_timeout_seconds=1,
    solve_timeout_seconds=1,
):
    if isinstance(expr, Eq):
        expr = expr.lhs - expr.rhs

    second_derivative = next(
        d for d in expr.atoms(Derivative)
        if len(d.variables) == 2
        and d.variables[0] == d.variables[1]
    )

    x = second_derivative.variables[0]
    r = regularize(expr)

    result = simpleKovacic(
        r,
        x,
        integration_timeout_seconds=integration_timeout_seconds,
        solve_timeout_seconds=solve_timeout_seconds,
    )
    return result

    


if __name__ == "__main__":
    x = symbols("x")
    y = Function("y")(x)

    print("\nRegularization examples")
    equations_to_regularize = [
        Eq(y.diff(x, 2) + 2*x*y.diff(x) + y, 0),
        3*y.diff(x, 2) + 6/x*y.diff(x) - 9*y,
    ]
    for equation in equations_to_regularize:
        r = regularize(equation)
        print(f"  input : {equation}")
        print(f"  r     : {r}")
        print(f"  output: z'' = ({r})*z")

    print("\nSecOrdSolver examples")
    equations_to_solve = [
        Eq(y.diff(x, 2) - y, 0),
        y.diff(x, 2) - y,
    ]
    for equation in equations_to_solve:
        result = SecOrdSolver(equation)
        print(f"  input   : {equation}")
        print(f"  solution: {result['solution']}")
