# -*- coding: utf-8 -*-
"""Prepare deterministic analysis inputs from committed fixed data.

Core analysis no longer depends on live Yahoo/yfinance downloads.  The four
monthly-return parts are the frozen 2016-09 to 2026-08 KRW-unhedged proxy
return sample originally sourced from Yahoo Finance.  This script concatenates
them, validates the committed correlation/covariance package, and writes the
single monthly_returns_krw.csv expected by downstream scripts.
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"

parts = sorted(DATA.glob("monthly_returns_krw_part*.csv"))
if len(parts) != 4:
    raise FileNotFoundError(f"Expected 4 frozen return parts, found {len(parts)}")

frames = []
for i, path in enumerate(parts):
    df = pd.read_csv(path)
    if i > 0:
        # Parts 2-4 intentionally omit the repeated header.
        if "Date" not in df.columns:
            df = pd.read_csv(path, header=None)
            df.columns = frames[0].columns
    frames.append(df)

monthly = pd.concat(frames, ignore_index=True)
monthly["Date"] = pd.to_datetime(monthly["Date"])
monthly = monthly.drop_duplicates("Date").sort_values("Date")
if len(monthly) != 120:
    raise ValueError(f"Expected 120 monthly observations, found {len(monthly)}")

proxy = pd.read_csv(DATA / "proxy_mapping.csv")
assets = proxy["asset"].tolist()
missing = [a for a in assets if a not in monthly.columns]
if missing:
    raise ValueError(f"Missing return columns: {missing}")

R = monthly[assets].astype(float)
rho = R.corr().values
sigma = proxy.set_index("asset").loc[assets, "sigma_cma"].astype(float).values
cov = np.outer(sigma, sigma) * rho

rho_ref = pd.read_csv(DATA / "correlation.csv", index_col=0).loc[assets, assets].astype(float).values
cov_ref = pd.read_csv(DATA / "covariance_cma.csv", index_col=0).loc[assets, assets].astype(float).values

corr_err = float(np.max(np.abs(rho - rho_ref)))
cov_err = float(np.max(np.abs(cov - cov_ref)))
if corr_err > 1e-12 or cov_err > 1e-12:
    raise ValueError(f"Frozen-data validation failed: corr_err={corr_err}, cov_err={cov_err}")

monthly.set_index("Date").to_csv(DATA / "monthly_returns_krw.csv", encoding="utf-8")
print(f"Prepared {len(monthly)} rows; corr_err={corr_err:.3e}; cov_err={cov_err:.3e}")
