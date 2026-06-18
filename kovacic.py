"""
Elementary function integration module - based on SymPy
Credits to https://github.com/sympy/sympy.git
"""

from sympy import *
from sympy.integrals import integrate
from sympy.integrals.risch import risch_integrate
from sympy.core.symbol import Symbol
from sympy.core.expr import Expr

__all__ = [
    'integrate',
    'risch_integrate',
    'integrate_elementary',
    'check_elementary',
    'Symbol',
]


def integrate_elementary(expr, variable, *args, **kwargs):
    """
    Calculate indefinite integral of elementary functions
    
    Parameters:
    -----------
    expr : sympy expression
    variable : sympy Symbol
    
    Returns:
    --------
    Indefinite integral result in sympy expression
    
    Examples:
    ---------
    >>> from sympy import symbols
    >>> x = symbols('x')
    >>> integrate_elementary(x**2, x)
    x**3/3
    
    >>> integrate_elementary(1/(1+x**2), x)
    atan(x)
    """
    return integrate(expr, variable, *args, **kwargs)


def check_elementary(expr, variable):
    """
    Check if an expression has an elementary integral
    
    Parameters:
    -----------
    expr : sympy expression
    variable : sympy Symbol
    
    Returns:
    --------
    bool
        Whether elementary integral exists
    """
    try:
        result = risch_integrate(expr, variable)
        return not isinstance(result, Integral)
    except Exception:
        return False


def integrate_trig(expr, variable):
    """Trigonometric function integration"""
    return integrate_elementary(expr, variable)


def integrate_exp(expr, variable):
    """Exponential function integration"""
    return integrate_elementary(expr, variable)


def integrate_rational(expr, variable):
    """Rational function integration"""
    return integrate_elementary(expr, variable)


def integrate_algebraic(expr, variable):
    """Algebraic function integration"""
    return integrate_elementary(expr, variable)


# Example usage
if __name__ == "__main__":
    x = symbols('x')
    
    print("=== SymPy Elementary Function Integration ===\n")
    
    # Example 1: Polynomial
    expr1 = x**2 + 2*x + 1
    print(f"∫({expr1})dx = {integrate_elementary(expr1, x)}\n")
    
    # Example 2: Trigonometric function
    expr2 = sin(x)
    print(f"∫(sin(x))dx = {integrate_elementary(expr2, x)}\n")
    
    # Example 3: Exponential function
    expr3 = exp(x)
    print(f"∫(e^x)dx = {integrate_elementary(expr3, x)}\n")
    
    # Example 4: Rational function
    expr4 = 1/(1 + x**2)
    print(f"∫(1/(1+x²))dx = {integrate_elementary(expr4, x)}\n")
    
    # Example 5: Check if an elementary integral exists
    expr5 = sin(x)/x
    print(f"sin(x)/x has elementary integral: {check_elementary(expr5, x)}\n")
