from sympy import *
from sympy import Poly, cancel, fraction, roots


'''
exactify(expr) converts input into an exact SymPy expression when possible.
For example, it turns decimal numbers into rationals.
For obvious numbers like 3.141592653589793, it will convert to pi.
'''
def exactify(expr):
    return nsimplify(sympify(expr))


'''
laurent_coeff(expr, x, c, k) computes the coefficient of (x-c)**k in the Laurent expansion of expr at c.
'''
def laurent_coeff(expr, x, c, k):
    expr = cancel(expr)
    y = Dummy("y")
    shifted = expr.subs(x, c + y)
    order = max(1, abs(k) + 3)
    ser = series(shifted, y, 0, order).removeO()

    return simplify(expand(ser).coeff(y, k))


'''
pole_coeff(expr, x, c, n) computes the coefficient of 1/(x-c)**n in the Laurent expansion of expr at c, which is the same as the coefficient of (x-c)**(-n).
'''
def pole_coeff(expr, x, c, n):
    return laurent_coeff(expr, x, c, -n)


'''
order_at_infinity(expr, x) computes the order of the pole at infinity.
'''
def order_at_infinity(expr, x):
    expr = cancel(expr)
    num, den = fraction(expr)
    if simplify(num) == 0:
        return oo
    
    num_poly = Poly(num, x)
    den_poly = Poly(den, x)

    return den_poly.degree() - num_poly.degree()


'''
get_poles_with_orders(r, x) returns finite poles and their corresponding orders.
This uses roots(den, x).
'''
def get_poles_with_orders(r, x):
    r = cancel(r)
    num, den = fraction(r)
    if simplify(num) == 0:
        return {}
    
    den_poly = Poly(den, x)
    poles = roots(den_poly.as_expr(), x)

    if sum(poles.values()) != den_poly.degree():
        raise NotImplementedError("Could not express all poles explicity using roots().")

    return poles


'''
givePoleAnalysis(r, x) analyzes the poles of the rational function r in C(x).
It inputs a rational function r and the variable x, and returns a dictionary containing:
- 'num': the numerator of r
- 'den': the denominator of r
- 'deg_num': degree of the numerator
- 'deg_den': degree of the denominator
- 'finite_poles': a dictionary of finite poles and their orders
- 'num_finite_poles': the number of finite poles
- 'pole_at_infinity': the order of the pole at infinity
- 'case_analysis': a dictionary indicating which cases of Kovacic algorithm are valid
'''
def givePoleAnalysis(r, x):
    r = exactify(r)

    try:
        r_cancelled = cancel(r)
        num, den = fraction(r_cancelled)

        # Special case: r = 0
        if simplify(num) == 0:
            return {
                "num": S.Zero,
                "den": S.One,
                "deg_num": -oo,
                "deg_den": S.Zero,
                "finite_poles": {},
                "num_finite_poles": 0,
                "pole_at_infinity": oo,
                "case_analysis": {
                    "case1_valid": True,
                    "case2_valid": False,
                    "case3_valid": False,
                },
            }

        num_poly = Poly(num, x)
        den_poly = Poly(den, x)

        deg_num = num_poly.degree()
        deg_den = den_poly.degree()
        O_inf = deg_den - deg_num

        finite_poles = get_poles_with_orders(r_cancelled, x)

    except Exception as e:
        raise ValueError(f"Error analyzing poles: {e}")

    # Case 1: No pole, or all finite poles have order 1 or even.
    #         If O(infinity) < 3, O(infinity) must be even.
    # Case 2: Has a finite pole of order 2 or odd order > 2.
    # Case 3: Has finite poles of order 1 or 2 only, and O(infinity) >= 2.
    # Case 4: None of the above.

    # Check Case 1
    case1_valid = True

    for pole_point, order in finite_poles.items():
        if order != 1 and order % 2 != 0:
            case1_valid = False
            break

    if case1_valid:
        if O_inf < 3 and O_inf % 2 != 0:
            case1_valid = False

    # Check Case 2: has pole of order 2 or odd order > 2
    case2_valid = False

    for pole_point, order in finite_poles.items():
        if order == 2 or (order > 2 and order % 2 != 0):
            case2_valid = True
            break

    # Check Case 3: has finite poles of order 1 or 2 only, O(infinity) >= 2
    case3_valid = True

    if len(finite_poles) == 0:
        case3_valid = False

    for pole_point, order in finite_poles.items():
        if order not in (1, 2):
            case3_valid = False
            break

    if O_inf < 2:
        case3_valid = False

    return {
        "num": num,
        "den": den,
        "deg_num": deg_num,
        "deg_den": deg_den,
        "finite_poles": finite_poles,
        "num_finite_poles": len(finite_poles),
        "pole_at_infinity": O_inf,
        "case_analysis": {
            "case1_valid": case1_valid,
            "case2_valid": case2_valid,
            "case3_valid": case3_valid,
        },
    }


'''
coeff_at_infinity(expr, x, k) computes the coefficient of x**k at infinity.
'''
def coeff_at_infinity(expr, x, k):
    expr = cancel(expr)
    y = Dummy("y")

    expr_y = expr.subs(x, 1/y)
    order = max(1, abs(k) + 5)

    ser = series(expr_y, y, 0, order).removeO()
    
    return simplify(expand(ser).coeff(y, -k))


'''
is_nonnegative_integer(expr) determines whether expr is a nonnegative integer.
Returns int(expr) if expr is a nonnegative integer.
Otherwise return None. 
'''
def is_nonnegative_integer(expr):
    expr = simplify(radsimp(expr))
    if expr.is_integer is True and expr.is_nonnegative is True:
        return int(expr)

    return None


'''
case3_E_order2_fromb(b, n) computes set E to pass to step 2.
'''
def case3_E_order2_fromb(b, n):
    root = sqrt(S.One + S(4)*b)
    E = []

    for k in range(-n//2, n//2 + 1):
        e = S(6) + S(12)*S(k)*root/S(n)
        e = simplify(radsimp(e))

        if e.is_integer is True:
            E.append(e)

    return sorted(set(E), key=lambda t: int(t))