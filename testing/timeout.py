# This file contains a robust timeout utility for long-running functions.
from concurrent.futures import ProcessPoolExecutor, TimeoutError
from typing import Callable, Any

def run_with_timeout(func: Callable, *args: Any, timeout_seconds: int) -> Any:
    """
    Runs a function in a separate process and raises a TimeoutError if it
    takes too long.

    Args:
        func: The function to execute.
        *args: The arguments to pass to the function.
        timeout_seconds: The timeout duration in seconds.

    Returns:
        The result of the function call.

    Raises:
        TimeoutError: If the function execution exceeds the timeout.
    """
    with ProcessPoolExecutor(max_workers=1) as executor:
        future = executor.submit(func, *args)
        try:
            return future.result(timeout=timeout_seconds)
        except TimeoutError as e:
            raise TimeoutError(f"Function {func.__name__} timed out after {timeout_seconds}s") from e
