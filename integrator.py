from multiprocessing import get_context

from sympy import solve
from sympy.integrals import integrate


def _integrate_worker(expr, x, queue):
    try:
        queue.put(("ok", integrate(expr, x)))
    except Exception as exc:
        queue.put(("error", repr(exc)))


def integrate_with_timeout(expr, x, timeout_seconds=10):
    context = get_context("spawn")
    queue = context.Queue()
    process = context.Process(target=_integrate_worker, args=(expr, x, queue))

    process.start()
    process.join(timeout_seconds)

    if process.is_alive():
        process.terminate()
        process.join()
        return None

    if queue.empty():
        return None

    status, payload = queue.get()
    if status == "error":
        raise RuntimeError(payload)

    return payload


def _solve_worker(equation, symbol, queue):
    try:
        queue.put(("ok", solve(equation, symbol)))
    except Exception as exc:
        queue.put(("error", repr(exc)))


def solve_with_timeout(equation, symbol, timeout_seconds=10):
    context = get_context("spawn")
    queue = context.Queue()
    process = context.Process(
        target=_solve_worker,
        args=(equation, symbol, queue),
    )

    process.start()
    process.join(timeout_seconds)

    if process.is_alive():
        process.terminate()
        process.join()
        return None

    if queue.empty():
        return None

    status, payload = queue.get()
    if status == "error":
        raise RuntimeError(payload)

    return payload
