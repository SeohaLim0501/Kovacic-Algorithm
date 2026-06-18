from sympy import *
from sympy.integrals import integrate
from sympy.integrals.risch import risch_integrate
from sympy.core.symbol import Symbol
from sympy.core.expr import Expr
from sympy import symbols, fraction, cancel, Poly, roots
from itertools import product


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

    # Check Case 3: has pole of order 1 or 2 only, O(infinity) >= 2
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
Kovacic's case 1, step 1 for finite poles. Outputs dictionaries:
- sqrtR
- alphaPlus
- alphaMinus
'''
def case1_finite_pole_data(r, x, finite_poles):
    r = exactify(r)
    r = cancel(r)

    sqrtR = {}
    alphaPlus = {}
    alphaMinus = {}
    bval = {}

    for c, order in finite_poles.items():
        if order == 1:
            sqrtR[c] = S.Zero
            alphaPlus[c] = S.One
            alphaMinus[c] = S.One

        elif order == 2:
            sqrtR[c] = S.Zero
            b = pole_coeff(r, x, c, 2)
            bval[c] = b
            alphaPlus[c] = simplify((S.One + sqrt(S.One + S(4)*b))/S(2))
            alphaMinus[c] = simplify((S.One - sqrt(S.One + S(4)*b))/S(2))

        elif order % 2 == 0:
            v = order // 2
            sqrt_expr = sqrt(r)
            sqrtR[c] = S.Zero
            coeffs = {}

            for i in range(2, v + 1):
                ai = laurent_coeff(sqrt_expr, x, c, -i)
                coeffs[i] = ai
                sqrtR[c] += ai / (x - c)**i

            sqrtR[c] = simplify(sqrtR[c])

            av = coeffs[v]
            if simplify(av) == 0:
                raise ValueError(f"Leading coefficient of sqrt(r) at pole {c} is zero.")
            
            b = simplify(pole_coeff(r, x, c, v + 1) - pole_coeff(sqrtR[c]**2, x, c, v + 1))
            bval[c] = b
            alphaPlus[c] = simplify((S(v) + b / av) / S(2))
            alphaMinus[c] = simplify((S(v) - b / av) / S(2))
        
        else: 
            raise ValueError(f"Unexpected pole order {order} at {c} in case 1 analysis.")

    return {
        "sqrtR": sqrtR,
        "alphaPlus": alphaPlus,
        "alphaMinus": alphaMinus,
    }


'''
Kovacic's case 1, step 1 for infinity pole. Outputs values:
- sqrtR_inf
- alphaPlus_inf
- alphaMinus_inf
'''
def case1_infinity_pole_data(r, x, O_inf):
    r = exactify(r)
    r = cancel(r)

    if O_inf > 2:
        return {
            "sqrtR_inf": S.Zero,
            "alphaPlus_inf": S.Zero,
            "alphaMinus_inf": S.One
        }
    
    if O_inf == 2:
        num, den = fraction(r)
        s_poly = Poly(num, x)
        t_poly = Poly(den, x)

        b = simplify(s_poly.LC() / t_poly.LC())
        root = sqrt(S.One + S(4) * b)
        
        return {
            "sqrtR_inf": S.Zero,
            "alphaPlus_inf": simplify((S.One + root) / S(2)),
            "alphaMinus_inf": simplify((S.One - root) / S(2))
        }
    
    if O_inf <= 0:
        if O_inf % 2 != 0:
            raise ValueError(f"Unexpected negative odd order at infinity")

        v = -O_inf // 2
        y = Dummy("y")
        sqrt_expr_y = sqrt(r.subs(x, 1/y))
        ser = series(sqrt_expr_y, y, 0, v + 3).removeO()
        sqrtR_inf = S.Zero
        coeffs = {}

        for i in range(0, v + 1):
            ai = expand(ser).coeff(y, -i)
            coeffs[i] = ai
            sqrtR_inf += ai * x**i
        
        sqrtR_inf = simplify(sqrtR_inf)
        av = coeffs[v]

        if simplify(av) == 0:
            raise ValueError("Leading coefficient of sqrt(r) at infinity is zero.")
        
        b1 = coeff_at_infinity(r, x, v - 1)
        b2 = coeff_at_infinity(sqrtR_inf**2, x, v - 1)
        b = simplify(b1 - b2)

        return {
            "sqrtR_inf": sqrtR_inf,
            "alphaPlus_inf": simplify((-S(v) + b / av) / S(2)),
            "alphaMinus_inf": simplify((-S(v) - b / av) / S(2))
        }
    
    raise ValueError("Unhandled case for infinity pole for case 1")
        
'''
Kovacic's case 1, step 2. Compute possible alpha[c] family to find omega.
Returns dictionary of:
- "d": d,
- "d_expr": d_expr,
- "omega": omega,
- "finite_signs": dict(zip(poles_list, finite_signs)),
- "inf_sign": inf_sign
'''
def case1_step2(r, x, finite_poles, finite_data, inf_data):
    sqrtR = finite_data['sqrtR']
    alphaPlus = finite_data['alphaPlus']
    alphaMinus = finite_data['alphaMinus']
    sqrtR_inf = inf_data['sqrtR_inf']
    alphaPlus_inf = inf_data['alphaPlus_inf']
    alphaMinus_inf = inf_data['alphaMinus_inf']

    poles_list = list(finite_poles.keys())
    candidates = []

    for finite_signs in product([1, -1], repeat=len(poles_list)):
        for inf_sign in [1, -1]:
            alpha_sum = S.Zero
            omega = S.Zero

            for c, sgn in zip(poles_list, finite_signs):
                if sgn == 1:
                    alpha_c = alphaPlus[c]
                else:
                    alpha_c = alphaMinus[c]
                
                alpha_sum += alpha_c
                omega += sgn * sqrtR[c] + alpha_c / (x - c)

            if inf_sign == 1:
                alpha_inf = alphaPlus_inf
            else:
                alpha_inf = alphaMinus_inf
            d_expr = simplify(radsimp(alpha_inf - alpha_sum))
            d = is_nonnegative_integer(d_expr)

            if d is None:
                continue

            omega += inf_sign * sqrtR_inf
            omega = simplify(radsimp(omega))

            candidates.append({
                "d": d,
                "d_expr": d_expr,
                "omega": omega,
                "finite_signs": dict(zip(poles_list, finite_signs)),
                "inf_sign": inf_sign
            })
            
    return candidates


'''
Kovacic's case 1, step 3. Compute possible polynomial p of degree d. Returns dictionary of:
- "p": p_sol,
- "omega": omega,
- "z": z,
- "z_evaluated": z_eval,
- "coeff_solution": sol
'''
def case1_step3(r, x, candidates):
    omega = simplify(candidates["omega"])
    r = exactify(r)
    r = cancel(r)
    d = candidates["d"]
    
    if d == 0:
        p = S.One
        coeff_symbols = []
    else:
        coeff_symbols = symbols(f"a0:{d}")
        p = x**d + sum(coeff_symbols[i] * x**i for i in range(d))

    expr = diff(p, x, 2) + 2*omega*diff(p, x) + (
        diff(omega, x) + omega**2 - r
    ) * p
    
    expr = cancel(together(expr))
    num, den = fraction(expr)
    num = expand(num)
    
    if simplify(num) == 0:
        z = simplify(p * exp(Integral(omega, x)))
        z_eval = simplify(p * exp(integrate(omega, x)))

        return {
            "p": p,
            "omega": omega,
            "z": z,
            "z_evaluated": z_eval,
            "coeff_soution": {}
        }

    poly = Poly(num, x, expand=True)
    equations = [Eq(c, 0) for c in poly.all_coeffs()]
    if not coeff_symbols:
        return None
    
    sol_list = solve(equations, coeff_symbols, dict=True)
    if not sol_list:
        return None
    
    for sol in sol_list:
        p_sol = simplify(p.subs(sol))
        check = diff(p_sol, x, 2) + 2*omega*diff(p_sol, x) + (
            diff(omega, x) + omega**2 - r
        ) * p_sol
        check = simplify(cancel(together(check)))
        if check == 0:
            z = simplify(p_sol * exp(Integral(omega, x)))
            Iomega = integrate(omega, x)
            z_eval = simplify(p_sol * exp(Iomega))

            return {
                "p": p_sol,
                "omega": omega,
                "z": z,
                "z_evaluated": z_eval,
                "coeff_solution": sol
            }

    return None


'''
simpleKovacic(r, x) implements Kovacic's algorithm for solving z'' = r*z.
It inputs a rational function r and the variable x, and returns the solution if it exists.
'''
def simpleKovacic(r, x):
    r = exactify(r)

    result = {
        'solution' : None,
        'status' : "Not implemented",
        'debug' : {}
    }

    if simplify(r) == 0:
        C1, C2 = symbols("C1 C2")
        result["solution"] = C1*x + C2
        result["status"] = "Solved trivial equation z'' = 0"
        return result
    
    pole_analysis = givePoleAnalysis(r, x)
    result['debug']['pole_analysis'] = pole_analysis

    finite_poles = pole_analysis['finite_poles']
    O_inf = pole_analysis['pole_at_infinity']

    if pole_analysis["case_analysis"]["case1_valid"]:
        finite_data = case1_finite_pole_data(r, x, finite_poles)
        inf_data = case1_infinity_pole_data(r, x, O_inf)
        candidates = case1_step2(r, x, finite_poles, finite_data, inf_data)
        
        result["debug"]["case1_finite_data"] = finite_data
        result["debug"]["case1_infinity_data"] = inf_data
        result["debug"]["case1_candidates"] = candidates

        if candidates:
            for cand in candidates:
                sol = case1_step3(r, x, cand)
                if sol is not None:
                    result["solution"] = sol["z_evaluated"]
                    result["status"] = "Solved by Kovacic case 1"
                    result["debug"]["successful_case"] = 1
                    result["debug"]["successful_candidate"] = cand
                    result["debug"]["case1_step3"] = sol
                    return result
        
        else:
            result["debug"]["case1_status"] = "Case 1 valid, but no nonnegative integer d found"

    else:
        result["debug"]["case1_status"] = "Case 1 conditions not satisfied"
        
    if pole_analysis["case_analysis"]["case2_valid"]:
        result["debug"]["case2_status"] = "Case 2 conditions satisfied, not implemented yet"

    if pole_analysis["case_analysis"]["case3_valid"]:
        result["debug"]["case3_status"] = "Case 3 conditions satisfied, not implemented yet"

    return result

if __name__ == "__main__":
    x = symbols('x')

    print("=" * 50)
    print("Simple Kovacic Algorithm")
    print("=" * 50 + "\n")

    test_cases = [
        # Keep these
        ("Test 1: z'' = -z", -1),
        ("Test 2: z'' = 0", 0),

        # Case 1 examples
        (
            "Test 3 [case 1]: z'' = (4*x**2 + 8*x + 6)/(2*x + 1)**2 * z",
            (4*x**2 + 8*x + 6)/(2*x + 1)**2
        ),
        (
            "Test 4 [case 1]: z'' = (7*x**2 + 10*x - 1)/(4*x**2*(x - 1)**4) * z",
            (7*x**2 + 10*x - 1)/(4*x**2*(x - 1)**4)
        ),
        (
            "Test 5 [case 1]: z'' = (3*x - 1)/(x**4*(x - 1)) * z",
            (3*x - 1)/(x**4*(x - 1))
        ),
        (
            "Test 6: z'' = (x**2 - 2*x + 3 + 1/x + 7/(4*x**2) - 5/x**3 + 1/x**4) * z",
            x**2 - 2*x + 3 + S.One/x + S(7)/(4*x**2) - S(5)/x**3 + S.One/x**4
        ),

        # Case 2 examples
        # These are shifted versions of the paper's case 2 example.
        # They have one order-2 pole and O(infinity)=1, so case 1 and case 3 fail.
        (
            "Test 7 [case 2]: z'' = (-8*x - 3)/(16*x**2) * z",
            (-8*x - 3)/(16*x**2)
        ),
        (
            "Test 8 [case 2]: z'' = (-8*x + 5)/(16*(x - 1)**2) * z",
            (-8*x + 5)/(16*(x - 1)**2)
        ),
        (
            "Test 9 [case 2]: z'' = (-8*x - 19)/(16*(x + 2)**2) * z",
            (-8*x - 19)/(16*(x + 2)**2)
        ),

        # Case 3 examples
        # Use these when testing the case 3 routine directly, or after disabling case 1 and case 2.
        # The first two come from the paper's worked case 3 section; the third is a shifted variant.
        (
            "Test 10 [case 3]: z'' = (4 - x)/(4*x*(x - 1)**2) * z",
            (4 - x)/(4*x*(x - 1)**2)
        ),
        (
            "Test 11 [case 3]: z'' = (24*x**2 + 40*x + 15)/(4*(x*(x + 1))**2) * z",
            (24*x**2 + 40*x + 15)/(4*(x*(x + 1))**2)
        ),
        (
            "Test 12 [case 3]: z'' = (24*x**2 - 8*x - 1)/(4*x**2*(x - 1)**2) * z",
            (24*x**2 - 8*x - 1)/(4*x**2*(x - 1)**2)
        ),
    ]

    for title, r in test_cases:
        print(title)
        result = simpleKovacic(r, x)
        print(result)
