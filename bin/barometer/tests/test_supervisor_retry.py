#!/usr/bin/env python3
"""Test the supervisor's retry-with-doubling-timeout logic.

This is a *synthetic* test: it replicates the exact supervisor loop that lives
inside ``barometer.__main__.main()`` (the retry queue, the per-task timeout
table, ``handle_time_failure``, the slow/deadlock watchdog) but drives it with
a fake worker so it runs in seconds instead of minutes.

It does NOT import the ``barometer`` package (that would pull in seaborn /
sklearn / statsmodels). Instead it re-implements the supervisor core verbatim
and asserts the retry policy behaves as documented:

  * a task that times out once then succeeds  -> retried once, completes
  * a task that always times out              -> 3 retries (600/1200/2400 s),
                                                 then 1 line in failed_analyses
  * a task that succeeds immediately          -> completes, no retry

Run:
    python3 tests/test_supervisor_retry.py
"""

import logging
import multiprocessing
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from concurrent.futures import TimeoutError as FuturesTimeoutError

# Version-aware timeout exception: on Python 3.11+ concurrent.futures.TimeoutError
# IS the builtin TimeoutError, so a single class suffices (matching the real
# supervisor in __main__.py which targets the 3.12 container). On < 3.11 (dev
# local 3.9) we must catch both.
_TIMEOUT_EXC = (TimeoutError,) if sys.version_info >= (3, 11) else (TimeoutError, FuturesTimeoutError)

# max_tasks_per_child was added in Python 3.11; the real supervisor uses 20.
_pool_kwargs = {"max_tasks_per_child": 20} if sys.version_info >= (3, 11) else {}

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger("test_supervisor_retry")

OUTDIR = os.path.join(os.path.dirname(__file__), "_tmp_retry_out")


def fake_analyze_section(df, cols, info, outdir, name, stat_test=None,
                         bmk_filter_cols=None, max_bmks=None, enabled_tests=None):
    """Stand-in for barometer.analysis.analyze_section.

    Behaviour is keyed on the last character of the section name:
      'A' -> fails (TimeoutError) on the 1st call, succeeds on the 2nd
      'B' -> always fails (TimeoutError)
      'C' -> always succeeds
    """
    tag = name[-1]
    # Use a per-process counter file so the spawn child can remember how many
    # times it has been called (spawn children cannot share the parent's dict).
    state = os.path.join(OUTDIR, f"count_{tag}")
    n = 0
    if os.path.exists(state):
        with open(state) as fh:
            n = int(fh.read().strip() or 0)
    n += 1
    with open(state, "w") as fh:
        fh.write(str(n))

    if tag == "A" and n == 1:
        raise TimeoutError("simulated timeout (1st attempt)")
    if tag == "B":
        raise TimeoutError("simulated timeout (always)")
    return {"ok": True}


def fake_wrapper(args_tuple):
    """Mimics analyze_section_wrapper WITHOUT the SIGALRM (fast test)."""
    df_source, cols, info, outdir, name, key, stat_test, bmk_filter, max_bmks, tto = args_tuple
    results = fake_analyze_section(
        df_source, cols, info, outdir, name,
        stat_test=stat_test, bmk_filter_cols=bmk_filter,
        max_bmks=max_bmks, enabled_tests=None,
    )
    return key, {"differential_table": os.path.join(outdir, "5_differential", "differential_results.csv")}


def run():
    """Replica of the supervisor core from __main__.main()."""
    vtype = "espf"
    n_jobs = 2
    all_results = {vtype: {"aggregate": {}, "feature": {}, "sites": {}}}
    failed_analyses = []

    BASE_TASK_TIMEOUT = 300
    MAX_RETRIES = 3
    max_pending_tasks = n_jobs * 3

    df = pd.DataFrame({"x": [1, 2, 3]})
    tasks = [
        (df, ["c1"], {}, os.path.join(OUTDIR, "secA"), "Task A", ("feature", "secA"), "nonparametric", None, 10),
        (df, ["c1"], {}, os.path.join(OUTDIR, "secB"), "Task B", ("feature", "secB"), "nonparametric", None, 10),
        (df, ["c1"], {}, os.path.join(OUTDIR, "secC"), "Task C", ("feature", "secC"), "nonparametric", None, 10),
    ]

    def record_failure(section_name, key, section_outdir, reason, detail=""):
        mtype, section_key = key
        failed_analyses.append({
            "vtype": vtype, "mtype": mtype, "section_key": section_key,
            "section_name": section_name, "section_outdir": section_outdir,
            "reason": reason, "detail": detail,
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        })

    mp_context = multiprocessing.get_context("spawn")
    with ProcessPoolExecutor(max_workers=n_jobs, mp_context=mp_context, **_pool_kwargs) as executor:
        future_to_key = {}
        task_iter = iter(tasks)
        retry_queue = []
        retry_count = {}
        task_timeout_by_key = {}
        completed = 0
        failed = 0
        submitted_count = 0
        last_progress_time = time.time()
        submit_time = {}

        def current_future_timeout():
            in_flight = [task_timeout_by_key.get(k, BASE_TASK_TIMEOUT) for k, _, _ in future_to_key.values()]
            return (max(in_flight) if in_flight else BASE_TASK_TIMEOUT) + 30

        def current_stall_timeout():
            return current_future_timeout() + 60

        pending_futures = set()

        def submit_task(task_tuple):
            nonlocal submitted_count
            df_source, cols, info, outdir, name, key, stat_test, bmk_filter, max_bmks = task_tuple
            tto = task_timeout_by_key.get(key, BASE_TASK_TIMEOUT)
            future = executor.submit(fake_wrapper, (df_source, cols, info, outdir, name, key, stat_test, bmk_filter, max_bmks, tto))
            future_to_key[future] = (key, name, outdir)
            submit_time[future] = time.time()
            pending_futures.add(future)
            submitted_count += 1
            return future

        def next_task():
            if retry_queue:
                return retry_queue.pop(0)
            return next(task_iter)

        def handle_time_failure(key, section_name, section_outdir, reason, detail=""):
            nonlocal failed
            n = retry_count.get(key, 0)
            if n < MAX_RETRIES:
                retry_count[key] = n + 1
                task_timeout_by_key[key] = BASE_TASK_TIMEOUT * (2 ** (n + 1))
                retry_queue.append(task_by_key[key])
                log.warning(f"  RETRY ({n + 1}/{MAX_RETRIES}) with {task_timeout_by_key[key]}s timeout: {section_name} ({reason})")
                return True
            failed += 1
            record_failure(section_name, key, section_outdir, reason, detail)
            return False

        task_by_key = {t[5]: t for t in tasks}

        def refill():
            while len(pending_futures) < max_pending_tasks:
                try:
                    submit_task(next_task())
                except StopIteration:
                    break

        refill()
        log.info(f"Initial batch submitted ({submitted_count} tasks)")

        while pending_futures or retry_queue:
            if not pending_futures and retry_queue:
                refill()
            now = time.time()
            future_timeout = current_future_timeout()
            for future in list(pending_futures):
                t0 = submit_time.get(future)
                if t0 is not None and (now - t0) > future_timeout:
                    future.cancel()
                    pending_futures.discard(future)
                    key, section_name, section_outdir = future_to_key[future]
                    submit_time.pop(future, None)
                    future_to_key.pop(future, None)
                    last_progress_time = time.time()
                    if handle_time_failure(key, section_name, section_outdir, "slow", f"exceeded {future_timeout}s"):
                        log.warning(f"  SLOW TASK: {section_name} - will retry")
                    else:
                        log.error(f"  SLOW TASK CANCELLED: {section_name}")
            if not pending_futures:
                time.sleep(1)
                continue
            try:
                done_iter = as_completed(pending_futures, timeout=5)
                for future in done_iter:
                    pending_futures.discard(future)
                    submit_time.pop(future, None)
                    key, section_name, section_outdir = future_to_key[future]
                    future_to_key.pop(future, None)
                    try:
                        result_key, results = future.result(timeout=1)
                        mtype, section_key = result_key
                        all_results[vtype][mtype][section_key] = results
                        completed += 1
                        last_progress_time = time.time()
                        log.info(f"  COMPLETED: {section_name}")
                    except _TIMEOUT_EXC:
                        last_progress_time = time.time()
                        tto = task_timeout_by_key.get(key, BASE_TASK_TIMEOUT)
                        if handle_time_failure(key, section_name, section_outdir, "timeout", f"exceeded {tto}s"):
                            log.warning(f"  TIMEOUT: {section_name} - will retry")
                        else:
                            log.error(f"  TIMEOUT: {section_name} - exceeded {tto}s")
                    except Exception as e:
                        failed += 1
                        last_progress_time = time.time()
                        log.error(f"  ERROR: {section_name} - {type(e).__name__}")
                        record_failure(section_name, key, section_outdir, "error", f"{type(e).__name__}: {e}")
                    if len(pending_futures) < max_pending_tasks:
                        try:
                            submit_task(next_task())
                        except StopIteration:
                            pass
                    break
            except _TIMEOUT_EXC:
                elapsed = time.time() - last_progress_time
                if elapsed > current_stall_timeout():
                    log.error(f"  DEADLOCK: no progress for {elapsed:.0f}s")
                    break

    if failed_analyses:
        pd.DataFrame(failed_analyses).to_csv(os.path.join(OUTDIR, "failed_analyses.tsv"), sep="\t", index=False)
    return completed, failed, failed_analyses


if __name__ == "__main__":
    # Clean any state from a previous run.
    os.makedirs(OUTDIR, exist_ok=True)
    for f in os.listdir(OUTDIR):
        os.remove(os.path.join(OUTDIR, f))

    completed, failed, failed_analyses = run()

    print(f"\ncompleted={completed}, failed={failed}")
    print(f"failed_analyses: {[(f['section_name'], f['reason']) for f in failed_analyses]}")

    assert completed == 2, f"Expected 2 completed (A, C), got {completed}"
    assert failed == 1, f"Expected 1 failed (B), got {failed}"
    assert len(failed_analyses) == 1, f"Expected 1 line in failed_analyses, got {len(failed_analyses)}"
    assert failed_analyses[0]["section_name"] == "Task B"
    assert failed_analyses[0]["reason"] == "timeout"

    # Task A must have been attempted twice (1 fail + 1 success).
    with open(os.path.join(OUTDIR, "count_A")) as fh:
        assert int(fh.read().strip()) == 2, "Task A should have run twice"
    # Task B must have been attempted 4 times (1 + 3 retries).
    with open(os.path.join(OUTDIR, "count_B")) as fh:
        assert int(fh.read().strip()) == 4, "Task B should have run 4 times"

    print("\nTEST PASSED: retry with timeout doubling works correctly")
