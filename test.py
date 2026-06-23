from argparse import ArgumentParser
from collections import defaultdict
from contextlib import redirect_stderr, redirect_stdout
import ctypes
import os
from pathlib import Path
import py_compile
import re
import sys
from time import perf_counter

from sympy import Derivative, Eq, Function, S, diff, symbols, sympify


TESTS_PATH = Path(__file__).with_name("kovacic_test.txt")
RESULTS_PATH = Path(__file__).with_name("result.txt")
PROJECT_FILES = (
    "kovacic.py",
    "utils.py",
    "integrator.py",
    "SecOrdSolver.py",
    "main.py",
)
GREEN = "\033[92m"
RED = "\033[91m"
RESET = "\033[0m"
ANSI_ESCAPE = re.compile(r"\x1b\[[0-9;]*m")


class TeeOutput:
    def __init__(self, console, result_file):
        self.console = console
        self.result_file = result_file

    def write(self, text):
        self.console.write(text)
        self.result_file.write(ANSI_ESCAPE.sub("", text))
        return len(text)

    def flush(self):
        self.console.flush()
        self.result_file.flush()


def enable_ansi_colors():
    if os.name != "nt":
        return

    try:
        kernel32 = ctypes.windll.kernel32
        stdout_handle = kernel32.GetStdHandle(-11)
        mode = ctypes.c_uint()
        if kernel32.GetConsoleMode(stdout_handle, ctypes.byref(mode)):
            kernel32.SetConsoleMode(stdout_handle, mode.value | 0x0004)
    except Exception:
        pass


def parse_args():
    parser = ArgumentParser(description="Validate Kovacic cases.")
    parser.add_argument(
        "-t",
        "--timeout",
        type=float,
        default=1,
        help="integration timeout in seconds (default: 1)",
    )
    parser.add_argument(
        "-t2",
        "--timer2",
        type=float,
        default=1,
        help="Case 3 omega solve timeout in seconds (default: 1)",
    )
    parser.add_argument(
        "-c",
        "--cases",
        type=Path,
        default=TESTS_PATH,
        help="path to the case test file (default: kovacic_test.txt)",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=RESULTS_PATH,
        help="path to the exported results (default: result.txt)",
    )
    parser.add_argument(
        "--start",
        type=int,
        default=1,
        help="1-based test number to start from (default: 1)",
    )
    return parser.parse_args()


def compile_project_files():
    project_dir = Path(__file__).parent
    for filename in PROJECT_FILES:
        py_compile.compile(
            str(project_dir / filename),
            doraise=True,
        )


def load_test_cases(path, x):
    test_cases = []
    y = Function("y")
    z = Function("z")
    locals_map = {
        "x": x,
        "y": y,
        "z": z,
        "S": S,
        "Eq": Eq,
        "Derivative": Derivative,
        "diff": diff,
    }

    text = path.read_text(encoding="utf-8-sig")
    for line_number, raw_line in enumerate(text.splitlines(), start=1):
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue

        parts = [part.strip() for part in line.split("|", 2)]
        if len(parts) != 3:
            raise ValueError(
                f"Line {line_number} must use: title | expression | expected_case"
            )

        title, expr_text, expected_text = parts
        expected_case = int(expected_text)
        if expected_case not in (0, 1, 2, 3, 4):
            raise ValueError(
                f"Line {line_number} has invalid expected case {expected_case}."
            )

        expression = sympify(expr_text, locals=locals_map)
        test_cases.append((title, expression, expected_case))

    return test_cases


def detected_case(result):
    if result.get("status") == "Solved trivial equation z'' = 0":
        return 0

    successful_case = result.get("debug", {}).get("successful_case")
    if successful_case in (1, 2, 3):
        return successful_case
    if result.get("status") == "No Liuvillian solution found":
        return 4
    return None


def print_mark(correct):
    if correct:
        print(f"Result: {GREEN}Correct!{RESET}")
    else:
        print(f"Result: {RED}Incorrect!{RESET}")


def run_tests(args):
    if args.timeout <= 0 or args.timer2 <= 0:
        raise ValueError("Timeout values must be positive.")
    if args.start <= 0:
        raise ValueError("--start must be a positive 1-based test number.")

    try:
        compile_project_files()
    except py_compile.PyCompileError as error:
        print(f"{RED}Incorrect!{RESET}")
        print(f"Compilation error: {error}")
        return 1

    from kovacic import simpleKovacic
    from main import print_test_result
    from SecOrdSolver import SecOrdSolver, regularize

    x = symbols("x")
    test_cases = load_test_cases(args.cases, x)
    if args.start > len(test_cases) + 1:
        raise ValueError(
            f"--start {args.start} is past the last test "
            f"({len(test_cases)})."
        )

    selected_cases = test_cases[args.start - 1:]
    correct_count = 0
    incorrect_count = 0
    incorrect_cases = []
    elapsed_by_case = defaultdict(list)

    print("=" * 72)
    print("Kovacic Case Validation")
    print(f"Input file: {args.cases}")
    print(f"Output file: {args.output}")
    print(f"Tests: {len(test_cases)}")
    print(f"Start test: {args.start}")
    print(f"Tests to run: {len(selected_cases)}")
    print(f"Integration timeout: {args.timeout:g} seconds")
    print(f"Equation solve timeout: {args.timer2:g} seconds")
    print("=" * 72)

    for index, (title, expr, expected_case) in enumerate(
        selected_cases,
        start=args.start,
    ):
        start = perf_counter()

        try:
            if expr.has(Derivative):
                r = regularize(expr)
                result = SecOrdSolver(
                    expr,
                    integration_timeout_seconds=args.timeout,
                    solve_timeout_seconds=args.timer2,
                )
            else:
                r = expr
                result = simpleKovacic(
                    r,
                    x,
                    integration_timeout_seconds=args.timeout,
                    solve_timeout_seconds=args.timer2,
                )

            elapsed_seconds = perf_counter() - start
            actual_case = detected_case(result)
            correct = actual_case == expected_case
            print_test_result(index, title, r, result, elapsed_seconds)
            print(f"Expected case: {expected_case}")
            print(f"Detected case: {actual_case if actual_case is not None else 'unknown'}")
            print_mark(correct)

        except Exception as error:
            elapsed_seconds = perf_counter() - start
            correct = False
            print("-" * 72)
            print(f"Test {index}: {title}")
            print(f"Elapsed: {elapsed_seconds:.3f} seconds")
            print(f"Expected case: {expected_case}")
            print(f"Error: {type(error).__name__}: {error}")
            print_mark(False)

        if correct:
            correct_count += 1
        else:
            incorrect_count += 1
            incorrect_cases.append(index)
        elapsed_by_case[expected_case].append(elapsed_seconds)

    total_count = len(selected_cases)
    score = (correct_count / total_count * 100) if total_count else 0.0
    print("=" * 72)
    print("Summary")
    print(f"Total: {total_count}")
    print(f"{GREEN}Correct: {correct_count}{RESET}")
    print(f"{RED}Incorrect: {incorrect_count}{RESET}")
    print(f"Score: {score:.2f}% ({correct_count}/{total_count})")
    if incorrect_cases:
        print(
            f"Incorrect test cases: {RED}"
            + ", ".join(str(index) for index in incorrect_cases)
            + RESET
        )
    else:
        print(f"Incorrect test cases: {GREEN}none{RESET}")
    print("Average elapsed by expected case:")
    for case_number in sorted(elapsed_by_case):
        case_times = elapsed_by_case[case_number]
        average_elapsed = sum(case_times) / len(case_times)
        case_label = "trivial" if case_number == 0 else f"Case {case_number}"
        print(
            f"  {case_label}: {average_elapsed:.3f} seconds "
            f"({len(case_times)} tests)"
        )
    print("=" * 72)

    return 0 if incorrect_count == 0 else 1


def main():
    enable_ansi_colors()
    args = parse_args()

    if args.cases.resolve() == args.output.resolve():
        raise ValueError("The input and output paths must be different.")

    with args.output.open("w", encoding="utf-8") as result_file:
        tee = TeeOutput(sys.stdout, result_file)
        with redirect_stdout(tee), redirect_stderr(tee):
            try:
                return run_tests(args)
            finally:
                from integrator import shutdown_symbolic_worker

                shutdown_symbolic_worker()


if __name__ == "__main__":
    raise SystemExit(main())
