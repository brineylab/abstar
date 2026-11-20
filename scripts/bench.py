#!/usr/bin/env python
"""Startup and throughput benchmark for the installed `abstar` CLI. Reports; never gates (exit 0)."""
import argparse, gzip, json, os, platform, resource, shutil, statistics, subprocess, sys, tempfile, time
from pathlib import Path

CORPUS = Path(__file__).resolve().parent.parent / "legacy/test_data/bcr_corpus/sequences.fasta.gz"
ABSTAR = str(Path(sys.executable).parent / "abstar")


def median_ms(args, n=10):
    ts = []
    for _ in range(n):
        t = time.perf_counter()
        subprocess.run([ABSTAR, *args], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        ts.append((time.perf_counter() - t) * 1000)
    return statistics.median(ts)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--replicate", type=int, default=20, help="corpus copies for the throughput run")
    ap.add_argument("--out", help="also write the JSON here")
    a = ap.parse_args()

    data = gzip.open(CORPUS, "rt").read()
    n_records = data.count(">") * a.replicate
    first_seq = data.split("\n", 2)[1]

    res = {"help_ms": median_ms(["--help"]), "single_seq_ms": median_ms(["run", "--seq", first_seq])}

    tmp = tempfile.mkdtemp(prefix="abstar-bench-")
    try:
        fa = Path(tmp) / "in.fasta"
        fa.write_text(data * a.replicate)
        out = Path(tmp) / "out"
        # Caveat: RUSAGE_CHILDREN ru_maxrss is the max over ALL children so far (incl. the
        # short runs above, which are smaller), so it reflects the throughput run here.
        t = time.perf_counter()
        subprocess.run([ABSTAR, "run", str(fa), str(out), "--quiet"], check=True)
        dt = time.perf_counter() - t
        res["throughput_rps"] = n_records / dt
        res["records"] = n_records
        # ru_maxrss is KiB on Linux, bytes on macOS.
        rss = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
        res["peak_rss_mb"] = rss / (1024 * 1024 if sys.platform == "darwin" else 1024)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    res["machine"] = {"platform": platform.platform(), "processor": platform.processor(),
                      "cpu_count": os.cpu_count(), "python": platform.python_version(),
                      "note": "not the spec's M-series laptop"}
    s = json.dumps(res, indent=2)
    print(s)
    if a.out:
        Path(a.out).write_text(s + "\n")
    ok_h, ok_s = res["help_ms"] <= 150, res["single_seq_ms"] <= 300
    print(f"help <=150ms: {'PASS' if ok_h else 'FAIL'} ({res['help_ms']:.0f}); "
          f"single_seq <=300ms: {'PASS' if ok_s else 'FAIL'} ({res['single_seq_ms']:.0f}); "
          f"throughput {res['throughput_rps']:.0f} rec/s (target >=10k on M-series, informational)")


main()
