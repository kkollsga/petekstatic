#!/usr/bin/env python3
"""Wait until a just-published package is installable from PyPI.

PyPI's JSON/project page can become visible before pip's simple index is usable
from every runner. Release workflows use this as the final PyPI visibility gate.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("package")
    parser.add_argument("version")
    parser.add_argument("--attempts", type=int, default=30)
    parser.add_argument("--sleep", type=float, default=20.0)
    parser.add_argument("--import-name", default=None)
    args = parser.parse_args()

    requirement = f"{args.package}=={args.version}"
    import_name = args.import_name or args.package.replace("-", "_")

    with tempfile.TemporaryDirectory(prefix="wait-pypi-") as tmp:
        for attempt in range(1, args.attempts + 1):
            target = Path(tmp) / f"site-{attempt}"
            target.mkdir()
            print(f"PyPI visibility check {attempt}/{args.attempts}: {requirement}", flush=True)
            install = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "pip",
                    "install",
                    "--disable-pip-version-check",
                    "--no-cache-dir",
                    "--only-binary=:all:",
                    "--target",
                    str(target),
                    requirement,
                ],
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
            )
            if install.returncode == 0:
                env = dict(os.environ)
                env["PYTHONPATH"] = str(target) + os.pathsep + env.get("PYTHONPATH", "")
                smoke = subprocess.run(
                    [sys.executable, "-c", f"import {import_name}; print({import_name}.__name__)"],
                    env=env,
                    text=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                )
                if smoke.returncode == 0:
                    print(smoke.stdout, end="")
                    return 0
                print(smoke.stdout, end="")
            else:
                print(install.stdout, end="")

            if attempt < args.attempts:
                time.sleep(args.sleep)

    print(f"{requirement} was not installable from PyPI in time", file=sys.stderr)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
