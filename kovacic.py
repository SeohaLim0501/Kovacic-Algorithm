"""
Elementary function integration module - based on SymPy
Credits to https://github.com/sympy/sympy.git
"""

from sympy import *
from sympy.integrals import integrate, indefinite_integrate
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
    주어진 함수가 초등함수 형태로 적분 가능한지 확인
    Check if an expression has an elementary integral
    
    Parameters:
    -----------
    expr : sympy expression
        검사할 함수 / Expression to check
    variable : sympy Symbol
        적분 변수 / Integration variable
    
    Returns:
    --------
    bool
        초등함수 적분 가능 여부 / Whether elementary integral exists
    """
    try:
        result = risch_integrate(expr, variable)
        # 결과가 Integral 객체가 아니면 초등함수로 표현 가능
        return not isinstance(result, Integral)
    except Exception:
        return False


# 간편한 사용을 위한 래퍼 함수들
def integrate_trig(expr, variable):
    """삼각함수 적분 / Trigonometric function integration"""
    return integrate_elementary(expr, variable)


def integrate_exp(expr, variable):
    """지수함수 적분 / Exponential function integration"""
    return integrate_elementary(expr, variable)


def integrate_rational(expr, variable):
    """유리함수 적분 / Rational function integration"""
    return integrate_elementary(expr, variable)


def integrate_algebraic(expr, variable):
    """대수함수 적분 / Algebraic function integration"""
    return integrate_elementary(expr, variable)


# 사용 예제
if __name__ == "__main__":
    x = symbols('x')
    
    print("=== SymPy Elementary Function Integration ===\n")
    
    # 예제 1: 다항식
    expr1 = x**2 + 2*x + 1
    print(f"∫({expr1})dx = {integrate_elementary(expr1, x)}\n")
    
    # 예제 2: 삼각함수
    expr2 = sin(x)
    print(f"∫(sin(x))dx = {integrate_elementary(expr2, x)}\n")
    
    # 예제 3: 지수함수
    expr3 = exp(x)
    print(f"∫(e^x)dx = {integrate_elementary(expr3, x)}\n")
    
    # 예제 4: 유리함수
    expr4 = 1/(1 + x**2)
    print(f"∫(1/(1+x²))dx = {integrate_elementary(expr4, x)}\n")
    
    # 예제 5: 초등함수 적분 가능 여부 확인
    expr5 = sin(x)/x
    print(f"sin(x)/x has elementary integral: {check_elementary(expr5, x)}\n")
