# -*- coding: utf-8 -*-
"""Wait until a marker appears in a log file, then run a command (used to chain studies without idle CPU)."""
import subprocess
import sys
import time

log, marker = sys.argv[1], sys.argv[2]
cmd = sys.argv[3:]
while True:
    try:
        if marker in open(log, encoding="utf-8", errors="ignore").read():
            break
    except FileNotFoundError:
        pass
    time.sleep(30)
sys.exit(subprocess.call(cmd))
