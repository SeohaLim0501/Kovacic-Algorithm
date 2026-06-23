from sympy import *
from sympy import Poly, cancel, fraction, roots


'''
exactify(expr) converts input into an exact SymPy expression when possible.
For example, it turns decimal numbers into rationals.
For obvious numbers like 3.141592653589793, it will convert to pi.
'''
def exactify(expr):
    expr = sympify(expr)
    if not expr.has(Float):
        return expr
    return nsimplify(expr)


'''
laurent_coeff(expr, x, c, k) computes the coefficient of (x-c)**k in the Laurent expansion of expr at c.
'''
def laurent_coeff(expr, x, c, k):
    expr = cancel(expr)
    y = Dummy("y")
    shifted = expr.subs(x, c + y)
    order = max(1, abs(k) + 3)

    try:
        ser = series(shifted, y, 0, order).removeO()
    except NotImplementedError:
        if not expr.is_rational_function(x):
            raise
        return _rational_laurent_coeff(expr, x, c, k)

    return normalize_complex_constant(expand(ser).coeff(y, k))


def _rational_laurent_coeff(expr, x, c, k):
    """Compute a rational Laurent coefficient using exact derivatives."""
    num, den = fraction(cancel(expr))
    pole_order = roots(Poly(den, x).as_expr(), x).get(c)

    if pole_order is None:
        if k < 0:
            return S.Zero
        return normalize_complex_constant(
            diff(expr, x, k).subs(x, c) / factorial(k)
        )

    coefficient_index = k + pole_order
    if coefficient_index < 0:
        return S.Zero

    q_leading = normalize_complex_constant(
        diff(den, x, pole_order).subs(x, c) / factorial(pole_order)
    )
    if q_leading == 0:
        raise ValueError(
            f"Could not determine the leading denominator coefficient at {c}."
        )

    coefficients = []
    for n in range(coefficient_index + 1):
        p_n = normalize_complex_constant(
            diff(num, x, n).subs(x, c) / factorial(n)
        )
        previous_terms = S.Zero

        for s in range(1, n + 1):
            q_term = normalize_complex_constant(
                diff(den, x, pole_order + s).subs(x, c)
                / factorial(pole_order + s)
            )
            previous_terms += q_term * coefficients[n - s]

        coefficients.append(
            normalize_complex_constant(
                (p_n - previous_terms) / q_leading
            )
        )

    return coefficients[coefficient_index]


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
    if num == 0:
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
    if num == 0:
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
        if num == 0:
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
    
    return normalize_complex_constant(expand(ser).coeff(y, -k))


'''
is_nonnegative_integer(expr) determines whether expr is a nonnegative integer.
Returns int(expr) if expr is a nonnegative integer.
Otherwise return None. 
'''
def is_nonnegative_integer(expr):
    expr = normalize_complex_constant(expr)
    if expr.is_integer is True and expr.is_nonnegative is True:
        return int(expr)

    return None


def normalize_complex_constant(expr):
    """Normalize an exact complex constant without touching symbolic expressions."""
    expr = cancel(sympify(expr))

    has_noninteger_power = any(
        power.exp.is_integer is not True
        for power in expr.atoms(Pow)
    )
    if not has_noninteger_power:
        return expr

    expr = radsimp(expr)

    if not expr.free_symbols:
        return simplify(expand_complex(expr))

    return expr


'''
case3_E_order2_fromb(b, n) computes set E to pass to step 2.
'''
def case3_E_order2_fromb(b, n):
    root = sqrt(S.One + S(4)*b)
    E = []

    for k in range(-n//2, n//2 + 1):
        e = S(6) + S(12)*S(k)*root/S(n)
        e = normalize_complex_constant(e)

        if e.is_integer is True:
            E.append(e)

    return sorted(set(E), key=lambda t: int(t))
