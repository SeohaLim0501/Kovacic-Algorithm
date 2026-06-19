from sympy import *
from sympy.integrals.risch import risch_integrate
from sympy.core.symbol import Symbol
from sympy.core.expr import Expr
from sympy import symbols, fraction, cancel, Poly
from itertools import product
from integrator import integrate_with_timeout, solve_with_timeout
from utils import (
    coeff_at_infinity,
    exactify,
    givePoleAnalysis,
    is_nonnegative_integer,
    laurent_coeff,
    pole_coeff,
    case3_E_order2_fromb,
)


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
- "inf_sign": inf_sign,
- "finite_alpha": selected alpha value for each finite pole,
- "infinite_alpha": selected alpha value at infinity
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
            finite_alpha = {}

            for c, sgn in zip(poles_list, finite_signs):
                if sgn == 1:
                    alpha_c = alphaPlus[c]
                else:
                    alpha_c = alphaMinus[c]

                finite_alpha[c] = alpha_c
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
                "inf_sign": inf_sign,
                "finite_alpha": finite_alpha,
                "infinite_alpha": alpha_inf
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
def case1_step3(r, x, candidates, integration_timeout_seconds=10):
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
        z = p * exp(Integral(omega, x))
        Iomega = integrate_with_timeout(omega, x, integration_timeout_seconds)
        integration_timed_out = Iomega is None
        if integration_timed_out:
            z_eval = z
        else:
            z_eval = simplify(p * exp(Iomega))

        return {
            "p": p,
            "omega": omega,
            "z": z,
            "z_evaluated": z_eval,
            "integration_timed_out": integration_timed_out,
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
            z = p_sol * exp(Integral(omega, x))
            Iomega = integrate_with_timeout(omega, x, integration_timeout_seconds)
            integration_timed_out = Iomega is None
            if integration_timed_out:
                z_eval = z
            else:
                z_eval = simplify(p_sol * exp(Iomega))

            return {
                "p": p_sol,
                "omega": omega,
                "z": z,
                "z_evaluated": z_eval,
                "integration_timed_out": integration_timed_out,
                "coeff_solution": sol
            }

    return None


'''
Kovasic's case 2, step 1 for finite poles. Returns dictionary of dictionaries:
- "E": E,
- "bval": bval
'''
def case2_finite_pole_data(r, x, finite_poles):
    r = exactify(r)
    r = cancel(r)
    E = {}
    bval = {}

    for c, order in finite_poles.items():
        if order == 1:
            E[c] = [S(4)]

        elif order == 2:
            b = pole_coeff(r, x, c, 2)
            bval[c] = b
            root = sqrt(S.One + S(4) * b)
            raw_E = [
                S(2),
                S(2) + S(2)*root,
                S(2) - S(2)*root
            ]

            Ec = []
            for e in raw_E:
                e = simplify(radsimp(e))
                if e.is_integer is True:
                    Ec.append(e)

            E[c] = sorted(set(Ec), key=lambda t: int(t))
        
        elif order > 2:
            E[c] = [S(order)]

        else:
            raise ValueError(f"Unexpected pole order {order} at {c} in case 2.")

    return {
        "E": E,
        "bval": bval
    }


'''
Kovasic's case 2, step 1 for infinite poles. Returns dictionary of dictionaries:
- "E_inf"
- "b_inf"
'''
def case2_infinite_pole_data(r, x, infinite_poles):
    r = exactify(r)
    r = cancel(r)

    O_inf = infinite_poles

    if O_inf > 2:
        return {
            "E_inf": [S(0), S(2), S(4)],
            "b_inf": None,
        }

    if O_inf == 2:
        num, den = fraction(r)
        num_poly = Poly(num, x)
        den_poly = Poly(den, x)

        b = simplify(num_poly.LC() / den_poly.LC())
        root = sqrt(S.One + S(4)*b)

        raw_E = [
            S(2),
            S(2) + S(2)*root,
            S(2) - S(2)*root,
        ]

        E_inf = []
        for e in raw_E:
            e = simplify(radsimp(e))
            if e.is_integer is True:
                E_inf.append(e)

        return {
            "E_inf": sorted(set(E_inf), key=lambda t: int(t)),
            "b_inf": b,
        }

    if O_inf < 2:
        return {
            "E_inf": [S(O_inf)],
            "b_inf": None,
        }

    raise ValueError("Unhandled infinity case in Kovacic case 2.")


'''
Kovasic's case 2, step 2, which computes possible d and corresponding theta. Returns list of dictionary of:
- "d"
- "d_expr"
- "theta"
- "e_finite"
- "e_inf"
'''
def case2_step2(r, x, finite_poles, finite_data, inf_data):
    E = finite_data["E"]
    E_inf = inf_data["E_inf"]
    poles_list = list(finite_poles.keys())
    candidates = []

    for c in poles_list:
        if len(E[c]) == 0:
            return []
        
    if len(E_inf) == 0:
        return []
    
    finite_choices = [E[c] for c in poles_list]

    for e_tuple in product(*finite_choices):
        for e_inf in E_inf:
            e_sum = S.Zero
            theta = S.Zero

            for c, e_c in zip(poles_list, e_tuple):
                e_sum += e_c
                theta += e_c / (x - c)

            theta = simplify(theta / S(2))

            d_expr = simplify((e_inf - e_sum) / S(2))
            d = is_nonnegative_integer(d_expr)

            if d is None:
                continue

            candidates.append({
                "d": d,
                "d_expr": d_expr,
                "theta": theta,
                "e_finite": dict(zip(poles_list, e_tuple)),
                "e_inf": e_inf
            })

    return candidates


'''
Kovasic's case 2, step 3. Computes possible polynomial p of degree d, expression phi, and omega to find the solution z. Returns dictionary:
"p": p_sol,
"theta": theta,
"phi": phi,
"omega": omega,
"z": z_formal,
"z_evaluated": z_eval,
"coeff_solution": sol,
"candidate": candidates
'''
def case2_step3(r, x, candidates, integration_timeout_seconds=10):
    theta = simplify(candidates["theta"])
    r = exactify(r)
    r = cancel(r)
    d = candidates["d"]
    
    if d == 0:
        p = S.One
        coeff_symbols = []
    else:
        coeff_symbols = symbols(f"a0:{d}")
        p = x**d + sum(coeff_symbols[i] * x**i for i in range(d))

    expr = (diff(p, x, 3) + 3*theta*diff(p, x, 2) + (3*theta**2 + 3*diff(theta, x, 1) - 4*r) * diff(p, x, 1)
    + (diff(theta, x, 2) + 3*theta*diff(theta, x, 1) + theta**3 - 4*r*theta - 2*diff(r, x, 1)) * p)
    
    expr = cancel(together(expr))
    num, den = expr.as_numer_denom()
    num = expand(num)
    
    if simplify(num) == 0:
        sol_list = [{}]

    else:
        poly = Poly(num, x, expand=True)
        equations = [Eq(c, 0) for c in poly.all_coeffs()]
        if not coeff_symbols:
            return None
        
        sol_list = solve(equations, coeff_symbols, dict=True)
        if not sol_list:
            return None
    
    for sol in sol_list:
        p_sol = simplify(p.subs(sol))
        check = (diff(p_sol, x, 3) + 3*theta*diff(p_sol, x, 2) + (3*theta**2 + 3*diff(theta, x, 1) - 4*r) * diff(p_sol, x, 1)
        + (diff(theta, x, 2) + 3*theta*diff(theta, x, 1) + theta**3 - 4*r*theta - 2*diff(r, x, 1)) * p_sol)
        check = simplify(cancel(together(check)))
        check_num, check_den = check.as_numer_denom()
        if simplify(expand(check_num)) != 0:
            continue
        
        phi = simplify(theta + diff(p_sol, x) / p_sol)
        constant_term = simplify((diff(phi, x) + phi**2 - 2*r) / S(2))
        discriminant = simplify(phi**2 - 4*constant_term)

        omega_candidates = [
            simplify((-phi + sqrt(discriminant)) / S(2)),
            simplify((-phi - sqrt(discriminant)) / S(2)),
        ]

        for omega in omega_candidates:
            z_formal = exp(Integral(omega, x))
            Iomega = integrate_with_timeout(omega, x, integration_timeout_seconds)
            integration_timed_out = Iomega is None
            if integration_timed_out:
                z_eval = z_formal
            else:
                z_eval = simplify(exp(Iomega))

            return {
                "p": p_sol,
                "theta": theta,
                "phi": phi,
                "omega": omega,
                "z": z_formal,
                "z_evaluated": z_eval,
                "integration_timed_out": integration_timed_out,
                "coeff_solution": sol,
                "candidate": candidates
            }

    return None


'''
Kovacic's case 3, step 1 for finite, infinite poles. Returns dictionary of:
- 4: {"E"/"E_inf":
      "bval"/"b_inf":}
- 6:
- 12:
'''
def case3_finite_pole_data(r, x, finite_poles):
    r = exactify(r)
    r = cancel(r)
    data = {}

    for n in [4, 6,12]:
        E = {}
        bval = {}

        for c, order in finite_poles.items():
            if order == 1:
                E[c] = [S(12)]
            
            elif order == 2:
                b = pole_coeff(r, x, c, 2)
                bval[c] = b
                E[c] = case3_E_order2_fromb(b, n)
            
            else:
                raise ValueError("This cannot happen in case 3")

        data[n] = {
            "E": E,
            "bval": bval
        }

    return data


def case3_infinite_pole_data(r, x, O_inf):
    r = exactify(r)
    r = cancel(r)
    if O_inf < 2:
        raise ValueError("Wrong allocation of case 3")
    
    if O_inf == 2:
        num, den = fraction(r)
        num_poly = Poly(num, x)
        den_poly = Poly(den, x)
        b = simplify(num_poly.LC() / den_poly.LC())

    else:
        b = S.Zero
    
    data = {}

    for n in [4, 6, 12]:
        E_inf = case3_E_order2_fromb(b, n)
        data[n] = {
            "E_inf": E_inf,
            "b_inf": b
        }

    return data


'''
Kovacic's case 3, step 2.
Returns a list of dictionaries:
- "n"
- "d"
- "d_expr"
- "theta"
- "S"
- "e_finite"
- "e_inf"
'''
def case3_step2(r, x, finite_poles, finite_data, inf_data):
    candidates = []
    poles_list = list(finite_poles.keys())

    S_poly = S.One
    for c in poles_list:
        S_poly *= (x - c)
    S_poly = expand(S_poly)

    for n in [4, 6, 12]:
        E = finite_data[n]["E"]
        E_inf = inf_data[n]["E_inf"]

        bad_n = False
        for c in poles_list:
            if len(E[c]) == 0:
                bad_n = True
                break

        if bad_n:
            continue

        if len(E_inf) == 0:
            continue

        finite_choices = [E[c] for c in poles_list]

        for e_tuple in product(*finite_choices):
            for e_inf in E_inf:
                e_sum = S.Zero
                theta = S.Zero

                for c, e_c in zip(poles_list, e_tuple):
                    e_sum += e_c
                    theta += e_c / (x - c)

                d_expr = simplify(S(n) * (e_inf - e_sum) / S(12))
                d = is_nonnegative_integer(d_expr)

                if d is None:
                    continue

                theta = simplify(S(n) * theta / S(12))

                candidates.append({
                    "n": n,
                    "d": d,
                    "d_expr": d_expr,
                    "theta": theta,
                    "S": S_poly,
                    "e_finite": dict(zip(poles_list, e_tuple)),
                    "e_inf": e_inf,
                })

    return candidates


'''
Kovacic's case 3, step 3. Compute p_{-1}.
'''
def case3_step3_find_P(r, x, candidate):
    r = exactify(r)
    r = cancel(r)
    n = candidate["n"]
    d = candidate["d"]
    theta = cancel(together(candidate["theta"]))
    S_poly = expand(candidate["S"])

    if d == 0:
        P = S.One
        coeff_symbols = []
    else:
        coeff_symbols = symbols(f"a0:{d}")
        P = x**d + sum(coeff_symbols[i] * x**i for i in range(d))

    S_theta = cancel(together(S_poly * theta))
    S2r = cancel(together(S_poly**2 * r))
    S_theta = expand(S_theta)
    S2r = expand(S2r)
    P_seq = {}
    P_seq[n] = -P
    P_seq[n + 1] = S.Zero

    for i in range(n, -1, -1):
        Pi = P_seq[i]
        Pi_plus = P_seq[i + 1]
        P_prev = (
            -S_poly * diff(Pi, x)
            + ((n - i) * diff(S_poly, x) - S_theta) * Pi
            - (n - i) * (i + 1) * S2r * Pi_plus
        )

        P_seq[i - 1] = cancel(together(P_prev))

    target = cancel(together(P_seq[-1]))
    target_num, target_den = target.as_numer_denom()
    target_num = expand(target_num)

    if simplify(target_num) == 0:
        P_sol = simplify(P)
        P_seq_sol = {
            k: simplify(v)
            for k, v in P_seq.items()
            if k <= n
        }

        return {
            "P": P_sol,
            "theta": theta,
            "S": S_poly,
            "n": n,
            "d": d,
            "P_sequence": P_seq_sol,
            "coeff_solution": {},
            "candidate": candidate
        }

    if not coeff_symbols:
        return None

    poly = Poly(target_num, x, expand=True)
    equations = [Eq(c, 0) for c in poly.all_coeffs()]

    sol_list = solve(
        equations,
        coeff_symbols,
        dict=True,
        simplify=False,
        rational=False
    )

    if not sol_list:
        return None
    
    for sol in sol_list:
        P_sol = simplify(P.subs(sol))
        target_check = simplify(target.subs(sol))
        target_check = cancel(together(target_check))
        check_num, check_den = target_check.as_numer_denom()

        if simplify(expand(check_num)) != 0:
            continue

        P_seq_sol = {
            k: simplify(v.subs(sol))
            for k, v in P_seq.items()
            if k <= n
        }

        return {
            "P": P_sol,
            "theta": theta,
            "S": S_poly,
            "n": n,
            "d": d,
            "P_sequence": P_seq_sol,
            "coeff_solution": sol,
            "candidate": candidate,
        }
    
    return None


'''
Kovacic's case 3, step 3. Compute polynomial for omega.
'''
def case3_step3_makeWeq(r, x, S_poly, P_sequence, n, omega_symbol = None):
    if omega_symbol is None:
        omega_symbol = Symbol("omega")
    
    S_poly = expand(S_poly)
    omega_eq = S.Zero
    
    for i in range(0, n + 1):
        Pi = P_sequence[i]
        omega_eq += S_poly**i * Pi * omega_symbol**i / factorial(n - i)

    omega_eq = cancel(together(omega_eq))
    omega_eq = collect(omega_eq, omega_symbol)

    poly_omega = Poly(omega_eq, omega_symbol, domain="EX")
    coeffs_by_degree = {
        degree_tuple[0]: coeff
        for degree_tuple, coeff in poly_omega.terms()
    }

    return {
        "omega_symbol": omega_symbol,
        "omega_equation": omega_eq,
        "omega_eq_zero": Eq(omega_eq, 0),
        "omega_poly": poly_omega,
        "omega_coeffs_by_degree": coeffs_by_degree,
    }

'''
simpleKovacic(r, x) implements Kovacic's algorithm for solving z'' = r*z.
It inputs a rational function r and the variable x, and returns the solution if it exists.
'''
def simpleKovacic(
    r,
    x,
    integration_timeout_seconds=10,
    solve_timeout_seconds=10,
):
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
                sol = case1_step3(r, x, cand, integration_timeout_seconds)
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
        finite_data = case2_finite_pole_data(r, x, finite_poles)
        inf_data = case2_infinite_pole_data(r, x, O_inf)
        candidates = case2_step2(r, x, finite_poles, finite_data, inf_data)
        
        result["debug"]["case2_finite_data"] = finite_data
        result["debug"]["case2_infinity_data"] = inf_data
        result["debug"]["case2_candidates"] = candidates

        if candidates: 
            for cand in candidates:
                sol = case2_step3(r, x, cand, integration_timeout_seconds)
                if sol is not None:
                    result["solution"] = sol["z_evaluated"]
                    result["status"] = "Solved by Kovacic case 2"
                    result["debug"]["successful_case"] = 2
                    result["debug"]["successful_candidate"] = cand
                    result["debug"]["case2_step3"] = sol
                    return result
                
        else: 
            result["debug"]["case2_status"] = "Case 2 valid, but no nonnegative integer d found"

    else:
        result["debug"]["case2_status"] = "Case 2 conditions not satisfied"

    if pole_analysis["case_analysis"]["case3_valid"]:
        finite_data = case3_finite_pole_data(r, x, finite_poles)
        inf_data = case3_infinite_pole_data(r, x, O_inf)
        candidates = case3_step2(r, x, finite_poles, finite_data, inf_data)
        
        result["debug"]["case3_finite_data"] = finite_data
        result["debug"]["case3_infinity_data"] = inf_data
        result["debug"]["case3_candidates"] = candidates

        if candidates:
            found_case3_P = False
            for cand in candidates:
                Presult = case3_step3_find_P(r, x, cand)
                if Presult is not None:
                    found_case3_P = True
                    result["debug"]["case3_N"] = cand["n"]
                    result["debug"]["case3_Efinite"] = cand["e_finite"]
                    result["debug"]["case3_Einfinite"] = cand["e_inf"]
                    result["debug"]["case3_P0"] = Presult["P_sequence"][0]
                    omega_eq_data = case3_step3_makeWeq(
                        r, x, Presult["S"], Presult["P_sequence"], Presult["n"],
                    )

                    omega_equation = omega_eq_data["omega_equation"]
                    omega_eq_zero = omega_eq_data["omega_eq_zero"]
                    omega_symbol = omega_eq_data["omega_symbol"]

                    omega_solutions = solve_with_timeout(
                        omega_eq_zero,
                        omega_symbol,
                        solve_timeout_seconds,
                    )
                    result["debug"]["case3_omega_equation"] = omega_equation

                    if omega_solutions is None:
                        result["debug"]["case3_solve_timed_out"] = True
                        result["debug"]["case3_status"] = (
                            "Case 3 omega solve timed out"
                        )
                    
                    else:
                        result["debug"]["case3_solve_timed_out"] = False
                        result["debug"]["case3_omega_solutions"] = omega_solutions
                        result["debug"]["case3_status"] = (
                            "Case 3 omega equation solved"
                        )

                    result["solution"] = {
                        "form": "exp(Integral(omega, x))",
                        "omega_symbol": omega_symbol,
                        "omega_equation": omega_eq_zero
                    }

                    result["status"] = f"Solved by Kovacic case 3 with n={cand['n']}"
                    result["debug"]["successful_case"] = 3
                    result["debug"]["successful_candidate"] = cand
                    result["debug"]["case3_step3"] = Presult

                    return result

            if not found_case3_P:
                result["debug"]["case3_status"] = (
                    "Case 3 valid, but no polynomial P found"
                )

        else:
            result["debug"]["case3_status"] = "Case 3 valid, but no nonnegative integer d found"

    else:
        result["solution"] = "No Liuvillian solution for this!"


    return result

