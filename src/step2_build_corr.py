# -*- coding: utf-8 -*-
"""STEP 2 - deterministic correlation/covariance build.

The assignment pipeline uses the committed frozen 120-month KRW proxy return
sample, not a live Yahoo/yfinance call.  Yahoo Finance remains the documented
original price source for the frozen sample, but network access is not required
for reproduction.

The full package build is implemented in team2_step2_corr_package.py.
"""
from pathlib import Path
import runpy

ROOT = Path(__file__).resolve().parents[1]
runpy.run_path(str(ROOT / "src" / "team2_step2_corr_package.py"), run_name="__main__")
