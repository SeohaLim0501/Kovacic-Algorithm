import atexit
from itertools import count
from multiprocessing import get_context
from queue import Empty
from threading import Lock

from sympy import solve
from sympy.integrals import integrate


_context = get_context("spawn")
_worker_process = None
_request_queue = None
_response_queue = None
_request_ids = count(1)
_worker_lock = Lock()


def _symbolic_worker(request_queue, response_queue):
    while True:
        request = request_queue.get()
        if request is None:
            return

        request_id, operation, args = request
        try:
            if operation == "integrate":
                payload = integrate(*args)
            elif operation == "solve":
                payload = solve(*args)
            else:
                raise ValueError(f"Unknown symbolic operation: {operation}")

            response_queue.put((request_id, "ok", payload))
        except Exception as exc:
            response_queue.put((request_id, "error", repr(exc)))


def start_symbolic_worker():
    global _worker_process, _request_queue, _response_queue

    with _worker_lock:
        if _worker_process is not None and _worker_process.is_alive():
            return

        _request_queue = _context.Queue()
        _response_queue = _context.Queue()
        _worker_process = _context.Process(
            target=_symbolic_worker,
            args=(_request_queue, _response_queue),
            daemon=True,
        )
        _worker_process.start()


def shutdown_symbolic_worker(force=False):
    global _worker_process, _request_queue, _response_queue

    with _worker_lock:
        process = _worker_process
        request_queue = _request_queue
        response_queue = _response_queue

        if process is None:
            return

        if process.is_alive() and not force:
            request_queue.put(None)
            process.join(timeout=2)

        if process.is_alive():
            process.terminate()
            process.join()

        for queue in (request_queue, response_queue):
            if queue is not None:
                queue.close()
                queue.join_thread()

        _worker_process = None
        _request_queue = None
        _response_queue = None


def _run_with_timeout(operation, args, timeout_seconds):
    start_symbolic_worker()

    with _worker_lock:
        request_id = next(_request_ids)
        _request_queue.put((request_id, operation, args))

        try:
            response_id, status, payload = _response_queue.get(
                timeout=timeout_seconds
            )
        except Empty:
            # Release the lock before shutdown_symbolic_worker acquires it.
            timed_out = True
        else:
            timed_out = False

        if not timed_out and response_id != request_id:
            raise RuntimeError("Received an unexpected symbolic worker response.")

    if timed_out:
        shutdown_symbolic_worker(force=True)
        return None

    if status == "error":
        raise RuntimeError(payload)

    return payload


def integrate_with_timeout(expr, x, timeout_seconds=1):
    return _run_with_timeout("integrate", (expr, x), timeout_seconds)


def solve_with_timeout(equation, symbol, timeout_seconds=1):
    return _run_with_timeout("solve", (equation, symbol), timeout_seconds)


atexit.register(shutdown_symbolic_worker)
