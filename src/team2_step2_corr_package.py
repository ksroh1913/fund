# -*- coding: utf-8 -*-
"""Rebuild Team 2 Step-2 correlation package from frozen source data.

Source provenance
-----------------
The frozen monthly-return sample was originally created from Yahoo Finance
prices for 069500.KS, URTH, EEM, 152380.KS, LQD, PSP, IGF, BIZD and KRW=X,
converted to KRW-unhedged returns for foreign assets.  The analytical pipeline
uses the committed frozen returns rather than a live API so a clean clone is
reproducible.  Live Yahoo regeneration is optional and not required to run the
assignment.

This script writes:
- data/monthly_returns_krw.csv
- data/correlation.csv
- data/covariance_cma.csv
- data/diagnostics.csv
- data/team2_step2_corr_package.xlsx
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
DATA.mkdir(exist_ok=True)

parts = sorted(DATA.glob("monthly_returns_krw_part*.csv"))
if len(parts) != 4:
    raise FileNotFoundError("Expected four committed monthly return parts.")

first = pd.read_csv(parts[0])
cols = first.columns.tolist()
frames = [first]
for p in parts[1:]:
    d = pd.read_csv(p, header=None)
    d.columns = cols
    frames.append(d)

monthly = pd.concat(frames, ignore_index=True)
monthly["Date"] = pd.to_datetime(monthly["Date"])
monthly = monthly.sort_values("Date").drop_duplicates("Date")
if len(monthly) != 120:
    raise ValueError(f"Expected 120 rows, got {len(monthly)}")

proxy = pd.read_csv(DATA / "proxy_mapping.csv")
assets = proxy["asset"].tolist()
R = monthly[assets].astype(float)
rho = R.corr()
sigma = proxy.set_index("asset").loc[assets, "sigma_cma"].astype(float)
cov = pd.DataFrame(
    np.outer(sigma.values, sigma.values) * rho.values,
    index=assets, columns=assets
)

corr_eig = np.linalg.eigvalsh(rho.values)
cov_eig = np.linalg.eigvalsh(cov.values)
diag = pd.DataFrame({
    "check": [
        "T_monthly","start_month","end_month",
        "corr_symmetry_max_abs_error","corr_diag_max_abs_error",
        "corr_min_eigenvalue","corr_condition_number",
        "cov_min_eigenvalue","cov_condition_number",
        "cov_diag_sigma_max_abs_error",
    ],
    "value": [
        len(monthly),
        monthly["Date"].min().date().isoformat(),
        monthly["Date"].max().date().isoformat(),
        float(np.max(np.abs(rho.values-rho.values.T))),
        float(np.max(np.abs(np.diag(rho.values)-1))),
        float(corr_eig.min()),
        float(corr_eig.max()/corr_eig.min()),
        float(cov_eig.min()),
        float(cov_eig.max()/cov_eig.min()),
        float(np.max(np.abs(np.sqrt(np.diag(cov.values))-sigma.values))),
    ]
})

monthly.set_index("Date").to_csv(DATA / "monthly_returns_krw.csv", encoding="utf-8")
rho.to_csv(DATA / "correlation.csv", encoding="utf-8")
cov.to_csv(DATA / "covariance_cma.csv", encoding="utf-8")
diag.to_csv(DATA / "diagnostics.csv", index=False, encoding="utf-8")

xlsx = DATA / "team2_step2_corr_package.xlsx"
with pd.ExcelWriter(xlsx, engine="openpyxl") as writer:
    proxy.to_excel(writer, sheet_name="proxy_mapping", index=False)
    monthly.set_index("Date").to_excel(writer, sheet_name="monthly_returns_krw")
    rho.to_excel(writer, sheet_name="correlation")
    cov.to_excel(writer, sheet_name="covariance_cma")
    diag.to_excel(writer, sheet_name="diagnostics", index=False)

print(f"Wrote {xlsx}")
