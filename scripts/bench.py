#!/usr/bin/env python
"""Startup and throughput benchmark for the installed `abstar` CLI. Reports; never gates (exit 0)."""
import argparse, gzip, json, os, platform, shutil, statistics, subprocess, sys, tempfile, time
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

    raw = gzip.open(CORPUS, "rb").read()
    n_records = raw.count(b">") * a.replicate
    first_seq = raw.split(b"\n", 2)[1].decode()

    res = {"help_ms": median_ms(["--help"]), "single_seq_ms": median_ms(["run", "--seq", first_seq])}

    tmp = tempfile.mkdtemp(prefix="abstar-bench-")
    try:
        fa = Path(tmp) / "in.fasta"
        with open(fa, "wb") as f:  # stream N copies; no large in-memory string
            for _ in range(a.replicate):
                f.write(raw)
        raw1 = raw
        # Run abstar under a tiny wrapper that reports its own RUSAGE_CHILDREN maxrss, so the
        # figure is abstar's alone and can't inherit this harness's high-water mark.
        wrap = ("import resource,subprocess,sys;subprocess.run(sys.argv[1:],check=True);"
                "print(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,file=sys.stderr)")
        # ru_maxrss is KiB on Linux, bytes on macOS.
        unit = 1024 * 1024 if sys.platform == "darwin" else 1024

        def run(path, n, threads=None):
            cmd = [sys.executable, "-c", wrap, ABSTAR, "run", str(path), str(Path(tmp) / "out"), "--quiet"]
            if threads:
                cmd += ["--threads", str(threads)]
            t = time.perf_counter()
            r = subprocess.run(cmd, check=True, stderr=subprocess.PIPE, text=True)
            return n / (time.perf_counter() - t), int(r.stderr.strip().splitlines()[-1]) / unit

        res["throughput_rps"], res["peak_rss_mb"] = run(fa, n_records)
        res["records"] = n_records
        # Single thread, and a one-copy input for the memory-flatness comparison (peak RSS vs input size).
        one = Path(tmp) / "one.fasta"
        one.write_bytes(raw1)
        res["throughput_1thread_rps"], res["peak_rss_1thread_mb"] = run(one, n_records // a.replicate, threads=1)
        res["peak_rss_1x_mb"] = run(one, n_records // a.replicate)[1]
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    res["machine"] = {"platform": platform.platform(), "processor": platform.processor(),
                      "cpu_count": os.cpu_count(), "python": platform.python_version(),
                      "note": "not the spec's M-series laptop"}
    s = json.dumps(res, indent=2)
    print(s)
    if a.out:
        Path(a.out).write_text(s + "\n")
    print(f"single thread {res['throughput_1thread_rps']:.0f} rec/s; peak RSS {res['peak_rss_1x_mb']:.0f} MB at {res['records'] // a.replicate} reads vs "
          f"{res['peak_rss_mb']:.0f} MB at {res['records']} reads")
    ok_h, ok_s = res["help_ms"] <= 150, res["single_seq_ms"] <= 300
    print(f"help <=150ms: {'PASS' if ok_h else 'FAIL'} ({res['help_ms']:.0f}); "
          f"single_seq <=300ms: {'PASS' if ok_s else 'FAIL'} ({res['single_seq_ms']:.0f}); "
          f"throughput {res['throughput_rps']:.0f} rec/s (target >=10k on M-series, informational)")


main()
