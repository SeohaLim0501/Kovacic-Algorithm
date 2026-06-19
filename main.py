from argparse import ArgumentParser
from pathlib import Path
from time import perf_counter

from sympy import S, sstr, symbols, sympify

from kovacic import simpleKovacic


EXAMPLES_PATH = Path(__file__).with_name("examples.txt")


def parse_args():
    parser = ArgumentParser(description="Run Kovacic algorithm examples.")
    parser.add_argument(
        "-t",
        "--timeout",
        type=float,
        default=10,
        help="integration timeout in seconds for case 1 and case 2 omega integration (default: 10)",
    )
    parser.add_argument(
        "-c",
        "--cases",
        type=Path,
        default=EXAMPLES_PATH,
        help="path to the examples file (default: examples.txt)",
    )
    parser.add_argument(
        "-t2",
        "--timer2",
        type=float,
        default=10,
        help="solve timeout in seconds for the case 3 omega equation (default: 10)",
    )
    return parser.parse_args()


def load_test_cases(path, x):
    test_cases = []
    locals_map = {"x": x, "S": S}

    for line_number, raw_line in enumerate(path.read_text().splitlines(), start=1):
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue

        if "|" in line:
            title, expr_text = [part.strip() for part in line.split("|", 1)]
        else:
            title = f"input line {line_number}"
            expr_text = line

        test_cases.append((title, sympify(expr_text, locals=locals_map)))

    return test_cases


def format_case_flags(case_analysis):
    labels = [
        ("case 1", case_analysis.get("case1_valid")),
        ("case 2", case_analysis.get("case2_valid")),
        ("case 3", case_analysis.get("case3_valid")),
    ]
    return ", ".join(f"{name}: {'yes' if valid else 'no'}" for name, valid in labels)


def format_finite_poles(finite_poles):
    if not finite_poles:
        return "none"

    return ", ".join(
        f"{sstr(pole)} (order {order})"
        for pole, order in finite_poles.items()
    )


def format_mapping(mapping):
    if not mapping:
        return "none"

    return ", ".join(
        f"{sstr(key)}: {sstr(value)}"
        for key, value in mapping.items()
    )


def format_expr_list(values):
    if not values:
        return "none"

    return "[" + ", ".join(sstr(value) for value in values) + "]"


def print_successful_case(debug):
    successful_case = debug.get("successful_case")
    successful_candidate = debug.get("successful_candidate")
    if not successful_candidate:
        return

    print(f"Case {successful_case} details:")

    if successful_case == 1:
        print("  alpha values:")
        print(
            "    finite_alpha: "
            f"{format_mapping(successful_candidate['finite_alpha'])}"
        )
        print(
            "    infinite_alpha: "
            f"{sstr(successful_candidate['infinite_alpha'])}"
        )

    elif successful_case == 2:
        print("  E sets:")
        finite_data = debug["case2_finite_data"]
        infinity_data = debug["case2_infinity_data"]
        for pole, values in finite_data["E"].items():
            print(f"    E_{sstr(pole)} = {format_expr_list(values)}")
        print(f"    E_infinity = {format_expr_list(infinity_data['E_inf'])}")

    elif successful_case == 3:
        print("  E sets:")
        n = successful_candidate["n"]
        finite_data = debug["case3_finite_data"][n]
        infinity_data = debug["case3_infinity_data"][n]
        for pole, values in finite_data["E"].items():
            print(f"    E_{sstr(pole)} = {format_expr_list(values)}")
        print(f"    E_infinity = {format_expr_list(infinity_data['E_inf'])}")

    print(f"  d: {successful_candidate['d']}")

    if successful_case == 1:
        print(f"  omega: {sstr(successful_candidate['omega'])}")

    elif successful_case == 2:
        step3 = debug.get("case2_step3", {})
        if step3:
            print(f"  omega: {sstr(step3['omega'])}")

    elif successful_case == 3:
        print(f"  N: {successful_candidate['n']}")
        print(f"  P_0: {sstr(debug['case3_P0'])}")


def print_test_result(index, title, r, result, elapsed_seconds):
    print("-" * 72)
    print(f"Test {index}: {title}")
    print(f"Elapsed: {elapsed_seconds:.3f} seconds")
    print(f"r(x): {sstr(r)}")
    print(f"Status: {result['status']}")

    solution = result.get("solution")
    if solution is not None:
        print(f"Solution: {sstr(solution)}")
    else:
        print("Solution: not found")

    debug = result.get("debug", {})
    pole_analysis = debug.get("pole_analysis")
    if pole_analysis:
        print("Pole analysis:")
        print(f"  finite poles: {format_finite_poles(pole_analysis['finite_poles'])}")
        print(f"  pole at infinity: order {pole_analysis['pole_at_infinity']}")
        print(f"Valid cases: {format_case_flags(pole_analysis['case_analysis'])}")

    print_successful_case(debug)


def main():
    args = parse_args()
    if args.timeout <= 0:
        raise ValueError("Timeout must be positive.")
    if args.timer2 <= 0:
        raise ValueError("Timeout2 must be positive.")

    x = symbols("x")
    examples_path = args.cases
    test_cases = load_test_cases(examples_path, x)

    print("=" * 72)
    print("Simple Kovacic Algorithm")
    print(f"Input file: {examples_path}")
    print(f"Integration timeout: {args.timeout:g} seconds")
    print(f"Equation Solve timeout: {args.timer2:g} seconds")
    print("=" * 72)

    for index, (title, r) in enumerate(test_cases, start=1):
        start = perf_counter()
        result = simpleKovacic(
            r,
            x,
            integration_timeout_seconds=args.timeout,
            solve_timeout_seconds=args.timer2,
        )
        elapsed_seconds = perf_counter() - start
        print_test_result(index, title, r, result, elapsed_seconds)


if __name__ == "__main__":
    main()
